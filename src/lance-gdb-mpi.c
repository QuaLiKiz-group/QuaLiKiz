
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

int lancegdbmpi_ (const char* argv0, int* lenargv0) {
    pid_t myPID;
    char cmyPID[1024];
    char argv0c[*lenargv0+1];
    strncpy(argv0c,argv0,*lenargv0);
    argv0c[*lenargv0+1]='\0';
    

    myPID=getpid();
    sprintf(cmyPID,"xterm -e gdb %s -p %d &",argv0c,myPID);
    printf( "Appel système : cmyPID=%s\n", cmyPID );
    system(cmyPID);
}

int attendgdbmpi_ () {
    system("sleep 2");
}
