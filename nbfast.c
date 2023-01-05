/* ---------------------------------------------------------------
Práctica 1-2.
Código fuente : nbfast.c
Grau Informàtica
48257747K   José Ramon Ariza Pérez
X9278777W   Mihaela Alexandra Buturuga
--------------------------------------------------------------- */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <stdio.h>
#include<unistd.h>
#include <string.h>
#include"pthread.h"
#include "stdbool.h"
#include "semaphore.h"

#ifdef D_GLFW_SUPPORT
#include<GLFW/glfw3.h>
#endif

// Macros to make code a little bit easier to understand because for speedup reasons, I'll use only 1D arrays
#define PX(i) (3*i+1)   //PosXY(i: Particle)
#define PY(i) (3*i+2)
#define MASS(i) (3*i+3)

#define VX(i) (4*i+0)   //VelocityXY(i: Particle)
#define VY(i) (4*i+1)
#define AX(i) (4*i+2)   //AccelerationXY(i: Particle)
#define AY(i) (4*i+3)

#define WIDTH 2000
#define HEIGHT 2000

#define BILLION 1000000000L

#define DEBUG 0

double G = 0.0001;
double dt = 0.005;
double rcutoff = 0.35;
double rlimit = 0.03;

//  Variables Globales
int nThread = 4;
int nShared = 500;      // Particles
int steps = 100;        // Iterations
int count;
int M = 25;

double *sharedBuff;     //Buffers to hold the position of the particles and their mass
double *localBuff;      //Buffer to hold velocity in x and y, and acceleration in x and y also
double *radius;
int *indexes;           //Index array, to speed up the creation of the tree (faster than passing the 3 floats per particle of x,y and mass)


int nLocal;
int nOriginal;

char filename[100];

struct Node *tree;

pthread_t *tid;

sem_t semTree, semIter;
pthread_barrier_t barr1;


bool graphicsEnabled = false;
bool inputFile = false;

void printPrueba(char *msg) {
    if (DEBUG) fprintf(stderr, "%s\n", msg);
}

struct Node {
    struct Node *children[4];   //Quaternary Tree
    int external;   // Determine if this node is external (far away enough)

    double CMX;     // Center Mass X,Y
    double CMY;
    double mass;
    double TRX;     // XY coordinate of top right corner
    double TRY;

    double LLX;     // XY coordinate of low left corner
    double LLY;

    double GCX;     // XY Node geometric center
    double GCY;
};

struct Statistics {
    double computTime;
    double totalComputTime;
    double loadImbalance;
    int evalPart;
    int deletedPart;
    int simplPart;
};

struct Partition {
    int id;
    int first;
    int last;
    struct Statistics stats;
};


struct Partition *attrs;
struct Statistics *st;

void freeMemory() {
    free(attrs);
    free(tid);
    free(sharedBuff);
    free(localBuff);
    free(indexes);
    free(radius);
}

void printError(const char *msg, int status) {
    fprintf(stderr, "%s\n", msg);
    freeMemory();
    exit(status);
}

void cancelRemainingThreads(int firstTid, int lastTid) {
    for (int i = firstTid; i < lastTid; i++) {
        fprintf(stderr, "Cancelling Thread %i\n", i);
        pthread_cancel(tid[i]);
    }
}

struct BuildTreeJob {
    struct Node *node;
    int *particleIndexes;
    int particleCount;
    int threads;
};

void initializeArgs(int argc, char **argv) {
    for (int j = 0; j < argc; ++j) {
        if (strcasecmp(argv[j], "-N") == 0)nShared = atoi(argv[j + 1]);
        else if (strcasecmp(argv[j], "-I") == 0) steps = atoi(argv[j + 1]);
        else if (strcasecmp(argv[j], "-T") == 0) nThread = atoi(argv[j + 1]);
        else if (strcasecmp(argv[j], "-M") == 0) M = atoi(argv[j + 1]);
        else if (strcasecmp(argv[j], "-G") == 0) { graphicsEnabled = true; }
        else if (strcasecmp(argv[j], "-F") == 0) {
            inputFile = true;
            strcpy(filename, argv[j + 1]);
        }
    }
}

void buildTreeConc(struct Node *node, int *indexes, int n, int numThreads);

void buildTreeThread(struct BuildTreeJob *job) {
    buildTreeConc(job->node, job->particleIndexes, job->particleCount, job->threads);
}

void buildTreeConc(struct Node *node, int *indexes, int n, int numThreads) {
    if (n == 1) { //This is an external node!
        node->external = 1;
        node->CMX = sharedBuff[PX(indexes[0])];
        node->CMY = sharedBuff[PY(indexes[0])];
        node->mass = sharedBuff[MASS(indexes[0])];
    } else {
        node->external = 0;
        //Arrays of indexes of particles per quartile
        int *NEi = (int *) malloc(sizeof(int) * n);
        int *NWi = (int *) malloc(sizeof(int) * n);
        int *SWi = (int *) malloc(sizeof(int) * n);
        int *SEi = (int *) malloc(sizeof(int) * n);
        int NWc = 0, SWc = 0, SEc = 0, NEc = 0;

        if (NEi == NULL || NWi == NULL || SWi == NULL || SEi == NULL) printError("Error Malloc", -1);

        int i;
        /** For each particle we will check where is it located relative to the geometric center,
            to sort them into the 4 children nodes**/
        for (i = 0; i < n; i++) {
            if (sharedBuff[PY(indexes[i])] < node->GCY) { //South half
                if (sharedBuff[PX(indexes[i])] < node->GCX) SWi[SWc++] = indexes[i];    //West wing
                else SEi[SEc++] = indexes[i];
            } else { //North half
                if (sharedBuff[PX(indexes[i])] < node->GCX) NWi[NWc++] = indexes[i];    //West wing
                else NEi[NEc++] = indexes[i];
            }
        }

        pthread_t *tidArray = 0;
        struct BuildTreeJob *jobArray = 0;
        int jobs, newThreads = 0, remainingThreads = 0, pendingThreads = 0, currentJob = 0;

        if (numThreads > 0) {   //Si la funcion tiene que crear threads, contar el trabajo a distribuir
            jobs = 0;
            if (NEc > 0) jobs++;
            if (NWc > 0) jobs++;
            if (SWc > 0) jobs++;
            if (SEc > 0) jobs++;

            newThreads = numThreads > jobs ? jobs
                                           : numThreads; //num de threads que se crearan en la funcion, como max 4, porque max 4 hijos
            remainingThreads = numThreads - newThreads; //num de threads que sobran, a repartir entre los hijos
            pendingThreads = newThreads;    //num de threads que quedan por crear en la funcion, al principio todos los newThreads
            currentJob = 0; //contador

            tidArray = malloc(sizeof(pthread_t) * newThreads);
            if (tidArray == NULL) printError("Error reservar vector de tids en buildTree", -2);
            jobArray = (struct BuildTreeJob *) malloc(sizeof(struct BuildTreeJob) * newThreads);
            if (jobArray == NULL) printError("Error reservar vector de jobs en buildTree", -3);
        }


        //If there are particles in the NorthWest quarter
        if (NEc > 0) {
            //This instruction declares a new node in the position 0
            node->children[0] = malloc(sizeof *node->children[0]);
            if (node->children[0] == NULL) printError("Error Malloc", -1);
            //We give the values of the Low Left and Top Right corner, and also the geometric center.
            node->children[0]->TRX = node->TRX;
            node->children[0]->TRY = node->TRY;
            node->children[0]->LLX = node->GCX;
            node->children[0]->LLY = node->GCY;
            node->children[0]->GCX = (node->GCX + node->TRX) / 2;
            node->children[0]->GCY = (node->GCY + node->TRY) / 2;

            //We build a tree in the new node, with the particles that are inside
            if (pendingThreads > 0) { //creando un nuevo thread si aun hay threads pendientes
                jobArray[currentJob].node = node->children[0];
                jobArray[currentJob].particleIndexes = NEi;
                jobArray[currentJob].particleCount = NEc;
                jobArray[currentJob].threads = remainingThreads / pendingThreads;

                if (pthread_create(&tidArray[currentJob], NULL, (void *(*)(void *)) buildTreeThread,
                                   (void *) &(jobArray[currentJob])) != 0) {
                    cancelRemainingThreads(0, i);
                    printError("Error crear hilo en buildTree", -4);
                }
                remainingThreads -= jobArray[currentJob].threads;
                pendingThreads--;
                currentJob++;

            } else buildTreeConc(node->children[0], NEi, NEc, 0);  //sino con una llamada recursiva
        } else node->children[0] = NULL;    //If not, we set the children to null

        //The next three blocks are exactly the same thing but for the other three nodes
        if (NWc > 0) {
            node->children[1] = malloc(sizeof *node->children[1]);
            if (node->children[1] == NULL) printError("Error Malloc", -1);
            node->children[1]->TRX = node->GCX;
            node->children[1]->TRY = node->TRY;
            node->children[1]->LLX = node->LLX;
            node->children[1]->LLY = node->GCY;
            node->children[1]->GCX = (node->LLX + node->GCX) / 2;
            node->children[1]->GCY = (node->GCY + node->TRY) / 2;

            if (pendingThreads > 0) { //creando un nuevo thread si aun hay threads pendientes
                jobArray[currentJob].node = node->children[1];
                jobArray[currentJob].particleIndexes = NWi;
                jobArray[currentJob].particleCount = NWc;
                jobArray[currentJob].threads = remainingThreads / pendingThreads;

                if (pthread_create(&tidArray[currentJob], NULL, (void *(*)(void *)) buildTreeThread,
                                   (void *) &(jobArray[currentJob])) != 0) {
                    cancelRemainingThreads(0, i);
                    perror("Error crear hilo en buildTree");
                }
                remainingThreads -= jobArray[currentJob].threads;
                pendingThreads--;
                currentJob++;
            } else buildTreeConc(node->children[1], NWi, NWc, 0);
        } else node->children[1] = NULL;

        if (SWc > 0) {
            node->children[2] = malloc(sizeof *node->children[2]);
            if (node->children[2] == NULL) printError("Error Malloc", -1);
            node->children[2]->TRX = node->GCX;
            node->children[2]->TRY = node->GCY;
            node->children[2]->LLX = node->LLX;
            node->children[2]->LLY = node->LLY;
            node->children[2]->GCX = (node->LLX + node->GCX) / 2;
            node->children[2]->GCY = (node->LLY + node->GCY) / 2;

            if (pendingThreads > 0) { //creando un nuevo thread si aun hay threads pendientes
                jobArray[currentJob].node = node->children[2];
                jobArray[currentJob].particleIndexes = SWi;
                jobArray[currentJob].particleCount = SWc;
                jobArray[currentJob].threads = remainingThreads / pendingThreads;

                if (pthread_create(&tidArray[currentJob], NULL, (void *(*)(void *)) buildTreeThread,
                                   (void *) &(jobArray[currentJob])) != 0) {
                    cancelRemainingThreads(0, i);
                    perror("Error crear hilo en buildTree");
                }
                remainingThreads -= jobArray[currentJob].threads;
                pendingThreads--;
                currentJob++;
            } else buildTreeConc(node->children[2], SWi, SWc, 0);
        } else node->children[2] = NULL;

        if (SEc > 0) {
            node->children[3] = malloc(sizeof *node->children[3]);
            if (node->children[3] == NULL) printError("Error Malloc", -1);
            node->children[3]->TRX = node->TRX;
            node->children[3]->TRY = node->GCY;
            node->children[3]->LLX = node->GCX;
            node->children[3]->LLY = node->LLY;
            node->children[3]->GCX = (node->GCX + node->TRX) / 2;
            node->children[3]->GCY = (node->LLY + node->GCY) / 2;

            if (pendingThreads > 0) { //creando un nuevo thread si aun hay threads pendientes
                jobArray[currentJob].node = node->children[3];
                jobArray[currentJob].particleIndexes = SEi;
                jobArray[currentJob].particleCount = SEc;
                jobArray[currentJob].threads = remainingThreads / pendingThreads;

                if (pthread_create(&tidArray[currentJob], NULL, (void *(*)(void *)) buildTreeThread,
                                   (void *) &(jobArray[currentJob])) != 0) {
                    cancelRemainingThreads(0, i);
                    perror("Error crear hilo en buildTree");
                }
                remainingThreads -= jobArray[currentJob].threads;
                pendingThreads--;
                currentJob++;

            } else buildTreeConc(node->children[3], SEi, SEc, 0);
        } else node->children[3] = NULL;

        node->mass = 0;
        node->CMX = 0;
        node->CMY = 0;

        //Esperar que los hijos acaben
        for (int c = 0; c < newThreads; c++)
            if (pthread_join(tidArray[c], NULL) != 0) {
                cancelRemainingThreads(c, newThreads);
                printError("Error Join", -1);
            }

        //Now that we have finished building the 4 trees beneath this node, we calculate the Center of Mass
        //based on the center of mass of the children
        for (i = 0; i < 4; i++) {
            if (node->children[i] != NULL) {
                node->mass += node->children[i]->mass;
                node->CMX += node->children[i]->CMX * node->children[i]->mass;
                node->CMY += node->children[i]->CMY * node->children[i]->mass;
            }
        }
        node->CMX = node->CMX / node->mass;
        node->CMY = node->CMY / node->mass;
        //And tadaaa
    }

}

void calculateForceNoParam(struct Node *_tree, int index, int *counter) {
    double distanceX = _tree->CMX - sharedBuff[PX(index)];
    double distanceY = _tree->CMY - sharedBuff[PY(index)];
    double distance = sqrt(distanceX * distanceX + distanceY * distanceY);
    //First we check if the node is not actually the same particle we are calculating
    if (distance > 0) {
        //Now, we know it is not because there is some distance between the Center of Mass and the particle
        //If the node is external (only contains one particle) or is far away enough, we calculate the force with the center of mass
        if (distance > rcutoff || _tree->external) {
            (*counter)++;
            double f;
            if (distance < rlimit) {
                f = G * _tree->mass / (rlimit * rlimit * distance);
            } else {
                f = G * _tree->mass / (distance * distance * distance);
            }
            localBuff[AX(index)] += f * distanceX;
            localBuff[AY(index)] += f * distanceY;
        } else {
            //If not, we recursively call the calculateForce() function in the children that are not empty.
            for (int s = 0; s < 4; s++) {
                if (_tree->children[s] != NULL) {
                    calculateForceNoParam(_tree->children[s], index, counter);
                }
            }
        }
    }
}

void moveParticleNoParam(struct Partition *p) {
    int first = p->first, last = p->last;
    for (int i = first; i < last; i++) {
        //Unprecise but fast euler method for solving the time differential equation
        double oldX = sharedBuff[PX(indexes[i])];
        double oldY = sharedBuff[PY(indexes[i])];
        sharedBuff[PX(indexes[i])] += localBuff[VX(indexes[i])] * dt + localBuff[AX(indexes[i])] * dt * dt * 0.5;
        sharedBuff[PY(indexes[i])] += localBuff[VY(indexes[i])] * dt + localBuff[AY(indexes[i])] * dt * dt * 0.5;
        localBuff[VX(indexes[i])] = (sharedBuff[PX(indexes[i])] - oldX) / dt;
        localBuff[VY(indexes[i])] = (sharedBuff[PY(indexes[i])] - oldY) / dt;
    }
}

#ifdef D_GLFW_SUPPORT
void drawParticle(double *shrdBuff, double *radius, int index){
    glBegin(GL_TRIANGLE_FAN);
    int k;
    glVertex2f(shrdBuff[PX(index)],shrdBuff[PY(index)]);
    for(k=0;k<20;k++){
        float angle=(float) (k)/19*2*3.141592;
        glVertex2f(shrdBuff[PX(index)]+radius[index]*cos(angle),shrdBuff[PY(index)]+radius[index]*sin(angle));
    }
    glEnd();
}

void drawBarnesHutDivisions(struct Node *rootNode){
    if(!rootNode->external){
        glBegin(GL_LINES);
        glVertex2f(rootNode->GCX,rootNode->LLY);
        glVertex2f(rootNode->GCX,rootNode->TRY);
        glVertex2f(rootNode->LLX,rootNode->GCY);
        glVertex2f(rootNode->TRX,rootNode->GCY);
        glEnd();
        int i;
        for(i=0;i<4;i++){
            if(rootNode->children[i]!=NULL){
                drawBarnesHutDivisions(rootNode->children[i]);
            }
        }
    }
}
#endif

void SaveGalaxy(int count, int nShared, int *indexes, double *sharedBuff);

void SaveGalaxyFile(char *filename, int nShared, int *indexes, double *sharedBuff);

void kickParticles();

void initializeRandomParticles();

void SaveGalaxy(int count, int _nShared, int *_indexes, double *_sharedBuff) {
    char _filename[100];
    sprintf(_filename, "./res/galaxy_%dB_%di.out", _nShared, count);
    SaveGalaxyFile(_filename, _nShared, _indexes, _sharedBuff);
}

void SaveGalaxyFile(char *_filename, int _nShared, int *_indexes, double *_sharedBuff) {
    int i;
    FILE *res = fopen(filename, "w");

    fprintf(res, "%d\n", _nShared);
    for (i = 0; i < _nShared; i++) {
        fprintf(res, "%d\t%e\t%e\t%e\n", _indexes[i], _sharedBuff[PX(_indexes[i])], _sharedBuff[PY(_indexes[i])],
                _sharedBuff[MASS(_indexes[i])]);
    }
    fclose(res);
}

void ReadGalaxyFileOptimized(char *_filename) {
    int ind;
    FILE *input;
    input = fopen(_filename, "r");
    if (input == NULL) printError("Error opening file.", 1);
    // Read number of bodies.
    if (fscanf(input, "%d\n", &nShared) < 1) printError("Error reading number of particles.", 1);

    // Reserve memory for indexes and particles.
    sharedBuff = (double *) malloc(sizeof(double) * (3 * nShared + 1));
    radius = (double *) malloc(sizeof(double) * (nShared));

    nLocal = nShared;
    nOriginal = nShared;

    indexes = (int *) malloc(sizeof(int) * nShared);
    localBuff = (double *) malloc(sizeof(double) * (4 * nShared));

    if (sharedBuff == NULL || radius == NULL || indexes == NULL || localBuff == NULL) printError("Error Malloc", -1);

    for (int j = 0; j < nShared; j++) {
        if (fscanf(input, "%d\t%le\t%le\t%le\n", &ind, &(sharedBuff[PX(j)]),
                   &(sharedBuff[PY(j)]), &(sharedBuff[MASS(j)])) < 4)
            printError("Error reading number of particles.", -1);
        indexes[j] = j;
        radius[j] = sqrt(sharedBuff[MASS(j)]) * 0.0025;

        localBuff[VX(j)] = 0;
        localBuff[VY(j)] = 0;
        localBuff[AX(j)] = 0;
        localBuff[AY(j)] = 0;
        //printf("Body %d: (%le,%le) %le\NBody_BarnesHut-Concurrente", ind,(*sharedBuff)[PX((*indexes)[j])],(*sharedBuff)[PY((*indexes)[j])],(*sharedBuff)[MASS((*indexes)[j])]);
    }
    fclose(input);
}

#define DSaveIntermediateState 1
#define DIntervalIntermediateState 100
#define DShowStatistics 1
#define DIntervalStatistics 1

clock_t StartTime, EndTime;
double TimeSpent;

void ShowWritePartialResults(int count, int _nOriginal, int _nShared, int *_indexes, double *_sharedBuff) {
//    if (DSaveIntermediateState && !(count % DIntervalIntermediateState))
//        SaveGalaxy(count, _nOriginal, _indexes, _sharedBuff);

    if (DShowStatistics && !(count % DIntervalStatistics)) {
        int j = 0;
        clock_t CurrentTime;
        CurrentTime = clock();
        TimeSpent = (double) (CurrentTime - StartTime) / CLOCKS_PER_SEC;
        //Mins = (int)TimeSpent/60;
        //Secs = (TimeSpent-(Mins*60));
        printf("[%.3f] Iteration %d => %d Bodies (%d) \t(Body %d: (%le, %le) %le).\n", TimeSpent, count, _nShared,
               _nOriginal, j, _sharedBuff[PX(_indexes[j])], _sharedBuff[PY(_indexes[j])],
               _sharedBuff[MASS(_indexes[j])]);
    }
}

void calculateForcesThread(struct Partition *attr) {
    int first = attr->first, last = attr->last;
    for (int j = first; j < last; j++) {
        //Set initial accelerations to zero
        localBuff[AX(indexes[j])] = 0;
        localBuff[AY(indexes[j])] = 0;
        // s = Node children index (Quaternary tree) Max: 4
        for (int s = 0; s < 4; s++) {
            //Recursively calculate accelerations
            if (tree->children[s] != NULL)
                calculateForceNoParam(tree->children[s], indexes[j], &attr->stats.simplPart);
        }
    }
}

void kickParticles() {
    for (int i = 0; i < nLocal; i++) {
        if (sharedBuff[PX(indexes[i])] <= 0 || sharedBuff[PX(indexes[i])] >= 1 || sharedBuff[PY(indexes[i])] <= 0 ||
            sharedBuff[PY(indexes[i])] >= 1) {
            nLocal--;
            for (int j = i; j < nLocal; j++) {
                indexes[j] = indexes[j + 1];
            }
            i--;
        }
    }
    nShared = nLocal;
}

void checkHelpMode(int argc, char **argv) {
    if (argc <= 1 || strcasecmp(argv[1], "-H") == 0) {
        printf("Usage: ./%s , -N <Number_of_Bodies> -i, -I <Number_of_iterations> -t, -T <Number_of_Threads> -M <Statistics_interval> [-f, -F <Path_to_Initial_Galaxy_File>] [-g, -G (set GUI Mode)]",
               argv[0]);
        exit(0);
    }
}

void initializeRandomParticles() {
    // Reserve memory for indexes and particles.
    sharedBuff = (double *) malloc(sizeof(double) * (3 * nShared + 1));
    radius = (double *) malloc(sizeof(double) * (nShared));

    nLocal = nShared;
    nOriginal = nShared;

    indexes = (int *) malloc(sizeof(int) * nShared);
    localBuff = (double *) malloc(sizeof(double) * (4 * nShared));

    if (sharedBuff == NULL || radius == NULL || indexes == NULL || localBuff == NULL) printError("Error Malloc", -1);

    for (int j = 0; j < nShared; j++) {
        // RANDOM POSITION AND MASS ARRAY DISTRIBUTION
        sharedBuff[PX(j)] = (float) (j) / (nShared - 1) * 0.8 + 0.1;
        sharedBuff[PY(j)] = (float) (rand() % 4096) / 4095 * 0.8 + 0.1;
        sharedBuff[MASS(j)] = (double) (rand() % 2048) / 2047 * 2 + 1;

        // GRAPHICAL RADIUS CALCULATION (DEPENDING ON MASS)
        radius[j] = sqrt(sharedBuff[MASS(j)]) * 0.0025;

        // PARTICLE NUMBER ARRAY
        indexes[j] = j;

        // SPEED AND ACCELERATION ARRAY (INIT TO 0)
        localBuff[VX(j)] = 0;
        localBuff[VY(j)] = 0;
        localBuff[AX(j)] = 0;
        localBuff[AY(j)] = 0;
    }
}

void printInputArgs(int part, int iter, int th, char *file) {
    printf("------------NBody With------------\n");
    printf("- Particles: %d\n- Iterations: %d\n- Threads: %d\n", part, iter, th);
    if (inputFile) printf("- File: %s\n", file);
    printf("----------------------------------\n");
}

double getElapsedTime(struct timespec startTime, struct timespec endTime) {
    double seconds = (double) (endTime.tv_sec - startTime.tv_sec);
    double nseconds = (double) (endTime.tv_nsec - startTime.tv_nsec);
    return seconds + nseconds / BILLION;
}

void printPartialStatistics(struct Partition *p) {
    struct Statistics s = p->stats;
    int particles = p->last - p->first;
    printf("Thread %i Statistics Iter %i ## ", p->id, count);
    printf("Part: %i [%i-%i]    Eval Part: %i   Deleted Part: %i    Mass Simpl: %i  ", particles, p->first, p->last - 1,
           s.evalPart, s.deletedPart, s.simplPart);
    printf("Thread Time: %.6f Total Time: %.6f  Average Iter Time: %.6f\n", s.computTime, s.totalComputTime,
           s.computTime / count);
}

void threadRoutine(struct Partition *attr) {
    printPrueba("Antes del WHILE en HIJO");
    struct Statistics *s = &attr->stats;

    while (count <= steps) {
        s->evalPart += attr->last - attr->first;
        printPrueba("Dentro del WHILE en HIJO");
        struct timespec startTime, endTime;
        if (clock_gettime(CLOCK_REALTIME, &startTime) < 0) printError("Error calculate startTime", -1);

        printPrueba("Antes del semáforo en HIJO");
        // Semáforo para esperar a la creación del árbol
        sem_wait(&semTree);

        printPrueba("Después del semáforo en HIJO");
        // Fuerzas
        calculateForcesThread(attr);
        // Mover partículas
        moveParticleNoParam(attr);
        // Barrera para hacer kick después
        int ret = pthread_barrier_wait(&barr1);
        if (ret != 0 && ret != PTHREAD_BARRIER_SERIAL_THREAD) {
            //Cancelar threads
        }

        printPrueba("Después de la barrera en HIJO");
        sem_wait(&semIter);
        if (clock_gettime(CLOCK_REALTIME, &endTime) < 0) printError("Error calculate endTime", -1);
        s->computTime = getElapsedTime(startTime, endTime);
        s->totalComputTime += s->computTime;
        if (count % M == 0) {
            printPartialStatistics(attr);
        }
    }
}

void sem_post_n(sem_t *sem, int n) {
    for (int i = 0; i < n; ++i) {
        sem_post(sem);
    }
}


#ifdef D_GLFW_SUPPORT
void graphicThread(GLFWwindow *window){
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glClear(GL_COLOR_BUFFER_BIT);

//            SaveGalaxy(count, nShared, indexes, sharedBuff);

    //This is only for visualization
    drawBarnesHutDivisions(tree);
    for(int k=0;k<nShared;k++){
        drawParticle(sharedBuff,radius,indexes[k]);
    }

//    t=glfwGetTime()-t;
//    if(t<0.013){
//        usleep(1000*1000*(0.013-t));
//    }

    glfwSwapBuffers(window);
    glfwPollEvents();
}
#endif


int main(int argc, char *argv[]) {
    checkHelpMode(argc, argv);
    StartTime = clock();

    // Overwrite default values if they are passed as arguments.
    initializeArgs(argc, argv);

    attrs = (struct Partition *) malloc(sizeof(struct Partition) * nThread);
    tid = (pthread_t *) malloc(sizeof(pthread_t) * (nThread));

    if (attrs == NULL || tid == NULL) printError("Error Malloc", -1);

    srand(time(NULL));

    if (inputFile && access(filename, F_OK) == 0) {
        /* Read bodies initial state from file */
        ReadGalaxyFileOptimized(filename);
    } else {
        initializeRandomParticles();
    }
    printInputArgs(nShared, steps, nThread, filename);

    //This is the main node, the one that holds the first four children nodes that make the calculation zone
    tree = malloc(sizeof *tree);
    if (tree == NULL) printError("Error Malloc", -1);

    //LLX is the x coordinate of the Low Left corner
    tree->LLX = 0;
    //This is the y coordinate..
    tree->LLY = 0;

    //Now the same but for the top right corner
    tree->TRX = 1;
    tree->TRY = 1;
    //The coordinates of the geometric center of the node in x and y
    tree->GCX = 0.5;
    tree->GCY = 0.5;

    // Save initial state.
    sprintf(filename, "./res/galaxy_%dB_initial.out", nOriginal);
    SaveGalaxyFile(filename, nShared, indexes, sharedBuff);

    count = 1;
    //If we need to visualize

#ifdef D_GLFW_SUPPORT
    if(graphicsEnabled){
        //If you only care about the algorithm, skip until next comment
        if(!glfwInit()){
            printf("Failed to start GLFW\NBody_BarnesHut-Concurrente");
            return -1;
        }
        GLFWwindow *window = glfwCreateWindow(WIDTH,HEIGHT,"Simulation",NULL,NULL);
        if(!window){
            printf("Failed to open window\NBody_BarnesHut-Concurrente");
            return -1;
        }
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0,1,0,1,0,1);
        glMatrixMode(GL_MODELVIEW);

        while(!glfwWindowShouldClose(window) && count<=steps){
            buildTreeConc(tree, indexes, nLocal, nThread);
            pthread_create(&tid[nThread-1], NULL, (void *(*)(void *)) graphicThread, window);

            if (nShared == nLocal && count > 1) {
                for (int j = 0; j < nThread - 1; j++) {
                    if (pthread_create(&tid[j], NULL, (void *(*)(void *)) calculateForcesThread, &attrs[j]) != 0) {
                        cancelRemainingThreads(0, j);
                        printError("ERROR CREATE", -3);
                    }
                }
                calculateForcesThread(&attrs[nThread - 1]);
            } else {
                nShared = nLocal;
                int firstNonAssigned = 0, currentJob = 0, remainingThreads = nThread, tasks = nLocal;
                for (int j = 0; j < nThread - 1; j++) {
                    currentJob = tasks / remainingThreads;
                    attrs[j].first = firstNonAssigned;
                    attrs[j].last = firstNonAssigned + currentJob;
                    if (pthread_create(&tid[j], NULL, (void *(*)(void *)) calculateForcesThread, &attrs[j]) != 0) {
                        cancelRemainingThreads(0, j);
                        printError("ERROR CREATE", -1);
                    }
                    firstNonAssigned += currentJob; //Maybe attrs[j].last
                    remainingThreads--;
                    tasks -= currentJob;
                }
                attrs[nThread - 1].first = firstNonAssigned;
                attrs[nThread - 1].last = firstNonAssigned + tasks;
                calculateForcesThread(&attrs[nThread - 1]);
            }
            for (int i = 0; i < nThread; i++) {
                if (pthread_join(tid[i], NULL) != 0) {
                    cancelRemainingThreads(i, nThread);
                    printError("ERROR JOIN", -2);
                }
            }

            // Move particles
            for (int j = 0; j < nThread - 1; j++) {
                if (pthread_create(&tid[j], NULL, (void *(*)(void *)) moveParticleNoParam, &attrs[j]) != 0) {
                    cancelRemainingThreads(0, j);
                    printError("ERROR CREATE", -3);
                }
            }
            moveParticleNoParam(&attrs[nThread-1]);
            for (int i = 0; i < nThread-1; i++) {
                if (pthread_join(tid[i], NULL) != 0) {
                    cancelRemainingThreads(i, nThread);
                    printError("ERROR JOIN", -4);
                }
            }

            // Kick out particle if it went out of the box (0,1)x(0,1)
            kickParticles();

            //To be able to store the positions of the particles
            ShowWritePartialResults(count, nOriginal, nLocal, indexes, sharedBuff);
            //We advance one step
            count++;
        }
        glfwTerminate();
    } else {
#endif
    //This is the pure algorithm, without visualization
    //system("mkdir res");

    if ((pthread_barrier_init(&barr1, NULL, nThread + 1)) < 0) printError("Error creating barrier", -1);

    sem_init(&semTree, 0, 0);
    sem_init(&semIter, 0, 0);

    count = 1;

    int firstNonAssigned = 0, currentJob = 0, remainingThreads = nThread, tasks = nLocal;
    //Asignación proporcional
    for (int j = 0; j < nThread; j++) {
        attrs[j].id = j;
        currentJob = tasks / remainingThreads;
        attrs[j].first = firstNonAssigned;
        attrs[j].last = firstNonAssigned + currentJob;
        if (pthread_create(&tid[j], NULL, (void *(*)(void *)) threadRoutine, &attrs[j]) != 0) {
            cancelRemainingThreads(0, j);
            printError("ERROR CREATE", -3);
        }
        firstNonAssigned += currentJob; //Maybe attrs[j].last
        remainingThreads--;
        tasks -= currentJob;
    }


    while (count <= steps) {
        nShared = nLocal;
        //First we build the tree
        buildTreeConc(tree, indexes, nLocal, nThread);

        firstNonAssigned = 0, currentJob = 0, remainingThreads = nThread, tasks = nLocal;
        for (int j = 0; j < nThread; j++) {
            attrs[j].id = j;
            currentJob = tasks / remainingThreads;
            attrs[j].first = firstNonAssigned;
            attrs[j].last = firstNonAssigned + currentJob;
            firstNonAssigned += currentJob; //Maybe attrs[j].last
            remainingThreads--;
            tasks -= currentJob;
        }
        sem_post_n(&semTree, nThread);
        printPrueba("Después del semáforo en MAIN");
        // Barrera
        int ret = pthread_barrier_wait(&barr1);
        if (ret != 0 && ret != PTHREAD_BARRIER_SERIAL_THREAD) {
            cancelRemainingThreads(0, nThread);
            //Cancelar threads
        }
        printPrueba("Después de la barrera en MAIN\n");

        // Kick out particle if it went out of the box (0,1)x(0,1)
        kickParticles();

        //To be able to store the positions of the particles
        ShowWritePartialResults(count, nOriginal, nLocal, indexes, sharedBuff);
        //We advance one step
        count++;
        sem_post_n(&semIter, nThread);
        printPrueba("Fin iteración en MAIN");
    }
#ifdef D_GLFW_SUPPORT
    }
#endif

    pthread_barrier_destroy(&barr1);
    sem_destroy(&semTree);
    sem_destroy(&semIter);
    EndTime = clock();
    TimeSpent = (double) (EndTime - StartTime) / CLOCKS_PER_SEC;
    printf("NBody Simulation took %.3f seconds.\n", TimeSpent);

    // Save final state.
    sprintf(filename, "./res/galaxy_%dB_%di_final.out", nOriginal, count - 1);
    SaveGalaxyFile(filename, nLocal, indexes, sharedBuff);
    freeMemory();
    return 0;
}


