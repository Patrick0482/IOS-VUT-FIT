#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <semaphore.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <limits.h>

// mapping for shared memory
#define MMAP(pointer) {(pointer) = mmap(NULL, sizeof(*(pointer)), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);}
#define UNMAP(pointer) {munmap( (pointer), sizeof((pointer)) );}


FILE *proj2out = NULL;

//semaphores

sem_t *sem_one = NULL;
sem_t *sem_oxygen = NULL;
sem_t *sem_hydrogen = NULL;
sem_t *sem_barrier_mutex = NULL;
sem_t *sem_turnstile1 = NULL;
sem_t *sem_turnstile2 = NULL;
sem_t *sem_line = NULL;


// shared variables

int *action_count = NULL;
int *oxygen_count = NULL;
int *hydrogen_count = NULL;
int *oxy_id = NULL;
int *hydro_id = NULL;
int *barrier_n = NULL;
int *molecul_count = NULL;
int NO, NH, TI, TB; 

//function init() used for initialization of everything

int init(){
    MMAP(action_count);
    MMAP(oxygen_count);
    MMAP(hydrogen_count);
    MMAP(oxy_id);
    MMAP(hydro_id);
    MMAP(barrier_n);
    MMAP(molecul_count);
 
    if ((sem_one = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;
    if ((sem_oxygen = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;
    if ((sem_hydrogen = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;
    if ((sem_barrier_mutex = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;
    if ((sem_turnstile1 = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;
    if ((sem_turnstile2 = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;
    if ((sem_line = mmap(NULL, sizeof(sem_t), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED, -1, 0)) == 0) return 1;

    sem_init(sem_one ,1 ,1);
    sem_init(sem_oxygen ,1 ,0);
    sem_init(sem_hydrogen, 1, 0);
    sem_init(sem_barrier_mutex, 1, 1);
    sem_init(sem_turnstile1, 1, 0);
    sem_init(sem_turnstile2, 1, 1);
    sem_init(sem_line ,1 ,1);

    return 0;
}

//function cleaner() used for cleaning the space mapped by the function init()

void cleaner(){
    UNMAP(action_count);
    UNMAP(oxygen_count);
    UNMAP(hydrogen_count);
    UNMAP(oxy_id);
    UNMAP(hydro_id);
    UNMAP(barrier_n);
    UNMAP(molecul_count);

    sem_destroy(sem_one);
    sem_destroy(sem_oxygen);
    sem_destroy(sem_hydrogen);
    sem_destroy(sem_barrier_mutex);
    sem_destroy(sem_turnstile1);
    sem_destroy(sem_turnstile2);
    sem_destroy(sem_line);

    munmap(sem_one,sizeof(sem_t));
    munmap(sem_oxygen,sizeof(sem_t));
    munmap(sem_hydrogen,sizeof(sem_t));
    munmap(sem_barrier_mutex,sizeof(sem_t));
    munmap(sem_turnstile1,sizeof(sem_t));
    munmap(sem_turnstile2,sizeof(sem_t));
    munmap(sem_line,sizeof(sem_t));

}

//function barrier() used for synchronization
//waiting for the third atom to go through and using the info that all 3 atoms went through
//also increasing molecul_count

void barrier(){
    sem_wait(sem_barrier_mutex);

        *barrier_n += 1;

        if ((*barrier_n) == 3){
            
            *molecul_count+=1;
            sem_wait(sem_turnstile2);
            sem_post(sem_turnstile1);
        }

    sem_post(sem_barrier_mutex);

    sem_wait(sem_turnstile1);
    sem_post(sem_turnstile1);

    sem_wait(sem_barrier_mutex);

        *barrier_n -= 1;

        if ((*barrier_n) == 0){
            sem_wait(sem_turnstile1);
            sem_post(sem_turnstile2);

        }

    sem_post(sem_barrier_mutex);

    sem_wait(sem_turnstile2);
    sem_post(sem_turnstile2);
}

//function oxygen() used for generating oxygens

void oxygen(int id){
    
    sem_wait(sem_one);
    *oxygen_count += 1;       
    *oxy_id += 1;

    sem_wait(sem_line);
    fprintf(proj2out,"%d: O %d: started\n",(*action_count)++,id);
    sem_post(sem_line);
   
    srand(time(0));                           //random time
    usleep(1000*(rand()%(TI+1)));
     
    sem_wait(sem_line);
    fprintf(proj2out,"%d: O %d: going to queue\n",(*action_count)++,id);
    sem_post(sem_line);
    // sem_post(sem_oxygen);

    if (NO*2 > NH && *oxy_id >= ((NH/2) + 1)) {             //checks if there is not enough H for O and writing it at the last O used
		sem_wait(sem_line);
			fprintf(proj2out, "%d: O %d: not enough H\n", (*action_count)++,id);
		sem_post(sem_line);

		sem_post(sem_one);
		return;
	}else{
        if(*hydrogen_count >=2){                //check for the amount we need for molecule 2-H, 1-O
            sem_post(sem_hydrogen);             //unlocking semaphore and letting through 2-H
            sem_post(sem_hydrogen);
            *hydrogen_count -= 2;

            sem_post(sem_oxygen);               //unlocking semaphore and letting through 1-O
            *oxygen_count -= 1;        
        }
        else{
            sem_post(sem_one);
        }
    }

    sem_wait(sem_oxygen);                   //waiting for the hydrogen thread to arrive

    barrier();                              //calling for the barrier function
    
    sem_wait(sem_line);
    fprintf(proj2out,"%d: O %d: creating molecule %d\n",(*action_count)++,id,*molecul_count);
    sem_post(sem_line);

    srand(time(0));                             
    usleep(1000*(rand()%(TB+1)));

    sem_wait(sem_line);
    fprintf(proj2out,"%d: O %d: molecule %d created\n",(*action_count)++,id,*molecul_count);
    sem_post(sem_line);
    
    sem_post(sem_one);
    
}

//function hydrogen() used for generating hydrogens

void hydrogen(int id){
    
    sem_wait(sem_one);
    *hydrogen_count += 1;
    *hydro_id += 1;

    sem_wait(sem_line);
    fprintf(proj2out,"%d: H %d: started\n",(*action_count)++,id);
    sem_post(sem_line);
     
    srand(time(0));
    usleep(1000*(rand()%(TI+1)));
     
    sem_wait(sem_line);
    fprintf(proj2out,"%d: H %d: going to queue\n",(*action_count)++,id);
    sem_post(sem_line);

    if (((NO*2) < NH && *hydro_id >= (NO*2)+1) || (NO*2 > NH && NH%2 && *hydro_id==NH)){    //checks if there is not enough H or O for H and writing it at the last H used
		sem_wait(sem_line);
			fprintf(proj2out, "%d: H %d: not enough O or H\n", (*action_count)++, id);
		sem_post(sem_line);

		sem_post(sem_one);
		return;
    }else{
        if(*hydrogen_count >=2 && *oxygen_count >= 1){     //similiar to the if in the oxygen() 
            sem_post(sem_hydrogen);                     //unlocking semaphore and letting through 2-H
            sem_post(sem_hydrogen);
            *hydrogen_count -= 2;

            sem_post(sem_oxygen);                       //unlocking semaphore and letting through 1-O
            *oxygen_count -= 1;        
        }else{
            sem_post(sem_one);
        }
    }
    
    sem_wait(sem_hydrogen);

    barrier();                      //calling for the barrier function

    sem_wait(sem_line);
    fprintf(proj2out,"%d: H %d: creating molecule %d\n",(*action_count)++,id,*molecul_count);
    sem_post(sem_line);

    srand(time(0));
    usleep(1000*(rand()%(TB+1)));

    sem_wait(sem_line);
    fprintf(proj2out,"%d: H %d: molecule %d created\n",(*action_count)++,id,*molecul_count);
    sem_post(sem_line);
    
}

int main(int argc, char *argv[])
{
    //shared variables from argv[]

    NO = atoi(argv[1]);
    NH = atoi(argv[2]);

    TI = atoi(argv[3]);
    if(TI < 0 || TI > 1000){                    //simple checks for bounds
        fprintf(stderr, "ERROR: TI argument out of bounds!\n");
        return 1;
    }
    
    TB = atoi(argv[4]);
    if(TB < 0 || TB > 1000){
        fprintf(stderr, "ERROR: TB argument out of bounds!\n");
        return 1;
    }

    //opening file / checking if it actually opened
    if ((proj2out = fopen("proj2.out","w")) == NULL){
        fprintf(stderr, "Error.\n");
        cleaner();
        exit(1);
    }
    
    //setting buffer of the output file to null for a organised output 
    setbuf(proj2out, NULL);

    //check for the argc amount
    if (argc != 5){
        fprintf(stderr, "ERROR: Wrong amount \n");
    }
    
    init();                 //initializing
    if(init() == 1){        //checks if initialization was successfull
        cleaner();
        fclose(proj2out);
        fprintf(stderr, "ERROR: Problem with semaphore inicialization!\n");
        exit(1);
    }

    //starting values for shared variables

    (*action_count) = 1;
    (*oxygen_count) = 0;
    (*hydrogen_count) = 0;
    (*oxy_id) = 0;
    (*hydro_id) = 0;
    (*barrier_n) = 0;
    (*molecul_count) = 0;

    //starting the oxygen thread

    for (int i = 1; i <= NO; i++){
         
        pid_t oxy = fork();
        if(oxy == 0) {
            oxygen(i);
            exit(0);
        }

        if(oxy < 0)
        {
            fprintf(stderr,"Fork failed!\n");       //checks if forking was successful
            exit(1);
        }
    }

    //starting the hydrogen thread 

    for (int i = 1; i <= NH; i++){
        pid_t hydro = fork();
        if(hydro == 0) {
            hydrogen(i);
            exit(0);
        }

        if(hydro < 0)
        {
            fprintf(stderr,"Fork failed!\n");       //checks if forking was successful
            exit(1);
        }
    }

    while(wait(NULL)>0);            //waiting for all semaphores to end            

    cleaner();                      //cleaning and closing the file
    fclose(proj2out);

    return 0;
}