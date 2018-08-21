#include <sys/types.h>
#if ( defined(Machine4) | defined(Machine5) | defined(MachineA) | defined(Machinel) | defined(Machinec) )

#endif

#if ( defined(Machinel) || defined(Machinep) )
#define TERMIOS
#endif


#include <string.h>


#ifdef access
#define close rmtclose
#define creat rmtcreat
#define dup rmtdup
#define fcntl rmtfcntl
#define fstat rmtfstat
#define mtio rmtmtio
#define isatty rmtisatty
#define lseek rmtlseek
#define lstat rmtlstat
#define open rmtopen
#define read rmtread
#define stat rmtstat
#define write rmtwrite

#ifndef Machinep
extern long rmtlseek ();	/* all the rest are int's */
#endif

#define setupoas rmtsetupoas
#define sendoas rmtsendoas
#define statf rmtstatf
#endif

#include <signal.h>
#include <sys/stat.h>
#include <sys/times.h>
#ifndef Machine5
#include <sys/ioctl.h>
#endif
#ifndef Machinel
#include <sys/mtio.h>
#else
#include <linux/mtio.h>
#endif
#include <sys/file.h>
#include <sys/socket.h>
#include <sys/errno.h>

#if ( defined(MachineT) ) 
#include <dirent.h>
#include <sys/fs/s5dir.h>
#define	_IO(x,y)	(IOC_VOID|('x'<<8)|y)
#define direct dirent
#define cgetenv_ CGETENV
#define clockf_ CLOCKF
#define cermes_ CERMES
#define cgterr_ CGTERR
#define csleep_ CSLEEP
#define ctrun_ CTRUN
#define cfork_ CFORK
#define cgetpid_ CGETPID
#define cgetppid_ CGETPPID
#define copendir_ COPENDIR
#define creaddir_ CREADDIR
#define cclosedir_ CCLOSEDIR
#define cmtio_ CMTIO
#define cread_ CREAD
#define clseek_ CLSEEK
#define copen_ COPEN
#define cunlink_ CUNLINK
#define cclose_ CCLOSE
#define cperror_ CPERROR
#define cwrite_ CWRITE
#define creadlink_ CREADLINK
#define cmkdir_ CMKDIR
#define cfstat_ CFSTAT
#define csetoas_ CSETOAS
#define cgethost_ CGETHOST
#define coassend_ COASSEND
#define coasnoretry_ COASNORETRY
#define usleep nousleep
#elif ( defined(MachineA) )
#include <dirent.h>
#include <s5fs/s5dir.h>
#define _IO(x,y)	(('x'<<8)|y)
#define direct dirent
#elif ( defined(Machine5) )
#include <dirent.h>
#define direct dirent
#elif ( defined(Machinec) )
#include <dirent.h>
#define direct dirent
#elif ( defined(Machinep) )
#include <dirent.h>
#define direct dirent
#else
#include <sys/dir.h>
#endif

#if ( defined(MachineT) )
#define F_STR char **
#define F_STRP(x) *x
#else
#define F_STR char *
#define F_STRP(x) x
#endif
#if ( defined(Machine5) | defined(Machinec) )
#include <sys/ttold.h>
#endif
#include <stdio.h>
#ifndef Machinep
#include <unistd.h>
#endif
#include <fcntl.h>
#ifndef Machinep
#include <sgtty.h>
#endif


#if ( defined(MachineS) | defined(MachineT) | defined(MachineA) )
#define	TIOCFLUSH	_IOW(t, 16, int)
#define	TIOCCDTR	_IO(t, 120)
#define	TIOCSDTR	_IO(t, 121)
#define	TIOCLGET	_IOR(t, 124, int)
#define	TIOCLSET	_IOW(t, 125, int)
#define	L_SET		0
#define	L_INCR		1
#define	L_XTND		2
#endif

#if ( defined(Machine5) | defined(Machinec)   )
#define	L_INCR		1
#define	L_XTND		2
#include <sys/filio.h>
#endif

#define LCBUF 100

typedef void Handler();
Handler *savehandler;
int oasichan;

#ifndef TERMIOS
static struct sgttyb   pack;
#else
#include <unistd.h>
#include <termios.h>
struct termios pack;
#endif

struct stat     buf;
static int oasintc;
static int timesup=0;
static int noretry=0;
DIR *dirpntr;

cabort_(str,lstr)
{       abort();
}

coasnoretry_()
{	noretry=1;
}



chektimesup_(value)
int *value;
{	*value = timesup;
}

exitfpe_()
{	extern void fpehandler();
	signal(SIGFPE,fpehandler);
}

void fpehandler(sig,code)
int sig, code;
{	printf("exiting on arithmetic exception\n");
	exit(2);
}

starttimesup_()
{	extern void settimesup();
	timesup=0;
 	signal(SIGALRM,settimesup);
}

void settimesup(sig,code)
int sig, code;
{	timesup=1;
}



csetintc_(value)
int *value;
{	oasintc = *value;
	printf("oasintc set to %d\n",oasintc);
}


#if ( defined(MachineT) | defined(MachineA) )
	cgethost_(mess,lmes)
#else
	cgethost_(mess,lmes,lmess)
	int lmess;
#endif
	F_STR mess;
	int *lmes;
{	extern int gethostname();
	int k, lmessi;
	char *namek;
#if ( defined(MachineT) | defined(MachineA) )
        F_STR temp;
        temp=mess ; lmessi=(int)(*(++temp));
#else
	lmessi=lmess;
#endif
	(void) gethostname(F_STRP(mess),lmessi);
	for ( k=0,namek=F_STRP(mess); ((*namek) != ' ') && ((*namek) != '\0') && (k < 80) ; k++,namek++)
	{	
	}
	*lmes=k;

}
     

#if ( defined(MachineT) )
	cgetenv_(name,mess,lmes)
#else
	cgetenv_(name,mess,lmes,lname,lmess)
	int lname, lmess;
#endif
	F_STR  name;
	F_STR  mess;
	int *lmes;
{	char buf[LCBUF];
 	char *i, *ptr;
	char *namek, *messk;
        int lnamei , lmessi;
	int k;
	extern char *getenv();
#if ( defined(MachineT) )
        F_STR temp;
        temp=name ; lnamei=(int)(*(++temp));
        temp=mess ; lmessi=(int)(*(++temp));
#else
	lnamei=lname;
	lmessi=lmess;
#endif
	for ( k=0,namek=F_STRP(name); ((*namek) != ' ') && ((*namek) != '\0') && k < lnamei ; )
	{	
		if(k >= LCBUF)
		{	printf("cgetenv: internal buffer size exceeded:%d\n",k);
			exit(2);
		}
		buf[k]=*namek; k++; namek++;
	}
	if(k >= LCBUF)
	{	printf("cgetenv: internal buffer size exceeded:%d\n",k);
		exit(2);
	}
	buf[k]='\0';
	ptr=getenv(buf);
	if( ptr == (char *)0 )
	{	*lmes=0;
	}
	else

	{
		for ( i=ptr, k=0, messk=F_STRP(mess); (*i) != '\0' ; )
		{	if ( k >= lmessi )
			{	printf("cgetenv: message buffer size exceeded:%d\n",k);
				exit(2);
			}
			(*messk)=(*i); i++; k++; messk++;
		}
		if ( k >= lmessi )
		{	printf("cgetenv: message buffer size exceeded:%d\n",k);
			exit(2);
		}
		(*messk)='\0';
		*lmes = k;
	}
}


clockf_(ichan,iopt,isize,ires,ierrno)
	int *ichan, *iopt, *isize, *ires, *ierrno;
{	
#ifndef Machinep
        extern int errno, lockf();
	int kopt;
        int tchan;
        int tsize, thou=1000;
	kopt = (int) (*iopt);
	
	switch (kopt) {
	case 0:
		kopt = F_ULOCK;
		break;
	case 1:
		kopt = F_LOCK;
		break;
	case 2:
		kopt = F_TLOCK;
		break;
	case 3:
		kopt = F_TEST;
		break;
	default:
		fprintf(stderr, "lockf: unknown option %d", *iopt);
		exit(2);
		break;
	}


	errno=0;
	tchan = (int)(*ichan);
	tsize = *isize;
/*	(void)cusleep_(&thou); */
/*	printf("calling lockf %d  %d  %d\n",tchan,kopt,tsize);  */
	*ires = (int) lockf(tchan, kopt, tsize);
/*	printf("lockf through\n");   */
/*	(void)cusleep_(&thou); */
	*ierrno = (int) errno;
#endif
}

cermes_(ierrno,mess,lmess,ll)
	char *mess;
	int *ierrno, *lmess, ll;
{       char *cp;
        cp=strerror(*ierrno);
        strcpy(mess,cp);
        *lmess=(int)strlen(mess);
}
	

copendir_(name,ichan)
	F_STR name;
	int *ichan;
{	extern DIR *opendir();
	dirpntr=opendir(F_STRP(name));
	*ichan=(int) dirpntr;
}

cclosedir_(ichan)
	int *ichan;
{
/*	closedir( (DIR *)(*ichan) ); */
	closedir( dirpntr );
}


creaddir_(ichan,name,lname,idno)
	F_STR name;
	int  *ichan, *lname, *idno;
{
	short i, l;
	struct direct *dptr;
	char *nn;
/*	if(dptr = readdir((DIR *)(*ichan)))     */
	if(dptr = readdir(dirpntr) )
	  {
/*	   l=dptr->d_namlen;
	   *lname= (int) l;
	   *idno=(int) dptr->d_ino;
           for (i=0, nn=name; i < l; i++, nn++)
             *nn=dptr->d_name[i];      old version */

	   *idno=(int) dptr->d_ino;
	   for (i=0, nn=F_STRP(name); dptr->d_name[i] != 0 ; i++, nn++)
		*nn=dptr->d_name[i];
	   *lname=i;
	   *(++nn)=dptr->d_name[i];
	  }
        else
          {
	   *lname = 0;
           *idno = -1;
          }
}

creadlink_(path, pbuf, pnbyt, ires, ierrno)
	F_STR pbuf;
        F_STR path;
	int           *pnbyt, *ires, *ierrno;
{
	extern int      errno;
	errno = 0;
	*ires = (int) readlink(F_STRP(path) , F_STRP(pbuf), (int) (*pnbyt));
	*ierrno = (int) errno;
}

cgterr_(ierrno)
	int *ierrno;
{	extern int errno;
	*ierrno=(int)errno;
}

cmkdir_(pname,mode,ires,ierrno)
        int *ires,*ierrno,*mode;
        F_STR pname;
{

        extern int errno, mkdir();
        errno=0;
        *ires=(int) mkdir(F_STRP(pname),(int)(*mode));
        *ierrno=(int)errno;
}
 


csetoas_(ichan,ierrno)
	int *ichan, *ierrno;
{	
	extern int setupoas();
	int jchan;
	jchan = *ichan;
	*ierrno = setupoas(jchan);
}

#if ( defined(MachineT) )
coassend_(ichan,scommand,sresponse,lr,iwt)
int *ichan, *lr, *iwt;
char **scommand, **sresponse;
{	char *command = *scommand, *response = *sresponse;
	int lencom = *((int *)(++scommand)), lenresp = *((int *)(++sresponse));
#else
coassend_(ichan,command,response,lr,iwt,lencom,lenresp)
int *ichan, *lr, *iwt, lencom, lenresp;
char *command, *response;
{
#endif
	extern int sendoas();
	(void)sendoas(*ichan,command,response,lr,*iwt,lencom,lenresp);
}


copen_(pname, ichan, iopt, ierrno, inew, mode)
	F_STR           pname;
	int           *ichan, *iopt, *ierrno, *inew, *mode;
{
	extern int      errno;
/*	extern int      open(); */
	int             kopt, knew;
	kopt = (int) (*iopt);
	knew = (int) (*inew);
	switch (kopt) {
	case 0:
		kopt = O_RDONLY;
		break;
	case 1:
		kopt = O_WRONLY;
		break;
	case 2:
		kopt = O_RDWR;
		break;
	case 8:
		kopt = O_APPEND;
		break;
	default:
		fprintf(stderr, "copen: unknown option %d", *iopt);
		exit(2);
		break;
	}
	switch (knew) {
	case 0:
		break;
	case 1:
		kopt = kopt | O_EXCL | O_CREAT;
		break;
	case 2:
		kopt = kopt | O_TRUNC | O_CREAT;
		break;
	case 3:
		kopt = kopt | O_TRUNC |  O_EXCL;
		break;
	default:
		fprintf(stderr, "copen: unknown status %d", *inew);
		exit(2);
		break;
	}
	errno = 0;
	*ichan = (int) open(F_STRP(pname), kopt|FNDELAY, (int)(*mode) );
	*ierrno = (int) errno;
}


cunlink_(pname,ires,ierrno)
	F_STR           pname;
	int *ires, *ierrno;
{
	extern int      errno;
	extern int      unlink();
	*ires = (int) unlink(F_STRP(pname));
	*ierrno = (int) errno;
}

cclose_(ichan, ires, ierrno)
	int           *ichan, *ires, *ierrno;
{
	extern int      errno;
	extern int      close();
	errno = 0;
	*ires = (int) close((int) (*ichan));
	*ierrno = (int) errno;
}

cread_(ichan, pbuf, pnbyt, ires, ierrno)
	char           *pbuf;
	int           *ichan, *pnbyt, *ires, *ierrno;
{
	extern int      errno;
	errno = 0;
	*ires = (int) read((int) (*ichan), pbuf, (int) (*pnbyt));
	*ierrno = (int) errno;
}

cwrite_(ichan, pbuf, pnbyt, ires, ierrno)
	char           *pbuf;
	int           *ichan, *pnbyt, *ires, *ierrno;
{
	extern int      errno;
	errno = 0;
	*ires = (int) write((int) (*ichan), pbuf, (int) (*pnbyt));
	*ierrno = (int) errno;
}

clseek_(ichan, offst, iopt, ires, ierrno)
	int           *ichan, *offst, *iopt, *ires, *ierrno;
{
	extern int      errno /*, lseek()  */;
#ifdef Machinep
	extern off_t lseek(int,off_t,int) ;
#else
	extern long lseek() ;
#endif
	int             kopt;
	kopt = (int) (*iopt);
	switch (kopt) {
	case 0:
		kopt = L_SET;
		break;
	case 1:
		kopt = L_INCR;
		break;
	case 2:
		kopt = L_XTND;
		break;
	default:
		fprintf(stderr, "clseek: unknown optiion %d", *iopt);
		exit(2);
		break;
	}
	errno = 0;
	*ires = (long) lseek((int) (*ichan), (*offst), kopt);
	*ierrno = (int) errno;
}


cmtio_(ichan, iop, icnt, ires, ierrno)
	int           *ichan, *iop, *icnt, *ires, *ierrno;
{
	extern int errno, mtio();
	errno=0;
	*ires = mtio(*ichan,*iop,*icnt);
	*ierrno = errno;
	if(errno != 0) (void)perror("cmtio");
}


cperror_(s)
	F_STR s;
{
	int             idum;
	perror(F_STRP(s));
}



csleep_(psec)
	int           *psec;
{
	int             idum;
	idum = sleep((unsigned) (*psec));
}

cfstat_(ichan, size, istat, ierrno)
	int           *ichan, *size, *istat, *ierrno;
{
	extern int errno;
	(void) statf(*ichan,size,istat);
	*ierrno = errno;
}

cftime_(path,ires,ierrno,itime,lpath)
	int *ires, *ierrno, *itime, *lpath ;
	char *path ;
{
	struct stat buf ;
	time_t modtim;
	extern int errno, stat() ;
	*ires = stat(path, &buf); 
	modtim = buf.st_mtime ; 
	*itime = modtim ; 
	*ierrno = errno ;
}

ctrun_(ichan, leng, ierrno)
	int           *ichan, *leng, *ierrno;
{
	extern int      errno, ftruncate();
	int             dummy;
	errno = 0;
	if ( *leng >= 0 ) {
	dummy = ftruncate((int) (*ichan), (unsigned long) (*leng));
	} else {
#if ( !defined(Machinel) && !defined(Machinep) )
	dummy = ftruncate((int) (*ichan), (unsigned long) tell(*ichan));
#else
	dummy = ftruncate((int) (*ichan), (unsigned long) ftell(*ichan));
#endif
	}
	*ierrno = (int) errno;
}

cfork_(ires,ierrno)
        int    *ires,*ierrno;
{
        extern int errno;
        errno = 0;
        *ires = (int) fork();
        *ierrno = (int) errno;
}

cgetpid_(ires,ierrno)
	int   *ires,*ierrno;
{
	extern int errno;
	errno = 0;
	*ires=(int)getpid();
	*ierrno = (int) errno;
}

cgetppid_(ires,ierrno)
	int   *ires,*ierrno;
{
	extern int errno;
	errno = 0;
	*ires=(int)getppid();
	*ierrno = (int) errno;
}


/* this used to be called creadlink - i don't think its used anywhere */
creadlinkg_(path, file, namlen, ierrno)
	char		*path, *file;
	int		*namlen, *ierrno;
{
	extern int	errno;
	int		dummy;
	errno = 0;
	dummy = readlink(path, file, namlen);
	if (dummy = -1) { *ierrno = errno; };
	if (dummy!= -1) { *ierrno = dummy; };
}

#if ( defined(Machine3) | defined(Machine4) | defined(Machine5) | defined(MachineA) | defined(Machinel) | defined(Machinep) | defined(Machinec) )
cutimes_(pname, times, ires, ierrno)
        int *ires, *ierrno;
        int *times;
        F_STR pname;
{
        extern int errno, utimes();
        errno=0;
        *ires= (int) utimes( F_STRP(pname),times);
        *ierrno=(int)errno;
}

cusleep_(psec)
	int           *psec;
{
	int             idum;
	idum = usleep((unsigned) (*psec));
}
#endif

#ifdef access

#undef access
#undef close
#undef creat
#undef dup
#undef fcntl
#undef fstat
#undef mtio
#undef isatty
#undef lseek
#undef lstat
#undef open
#undef read
#undef stat
#undef write

#undef setupoas
#undef sendoas
#undef statf

#endif
#if ( defined(Machine4) | defined(Machine5) | defined(MachineT) | defined(MachineA) | defined(Machinel) | defined(Machinep) | defined(Machinec))
int mtio(ichan, iop, icnt)
	int ichan, iop, icnt;
{
	extern int      errno;
	int             kopt, rc;
	struct mtop     magop;
	struct mtget     magstat;
	switch (iop) {
	case 0:
		kopt = MTWEOF;
		break;
	case 1:
		kopt = MTFSF;
		break;
	case 2:
		kopt = MTBSF;
		break;
	case 3:
		kopt = MTFSR;
		break;
	case 4:
		kopt = MTBSR;
		break;
	case 5:
		kopt = MTREW;
		break;
	case 6:
		kopt = MTOFFL;
		break;
	case 7:
		kopt = MTNOP;
		break;
	default:
		fprintf(stderr, "cmtio: unknown optiion %d", iop);
		exit(2);
		break;
	}
	errno = 0;
	magop.mt_op = kopt;
	magop.mt_count = icnt;
	rc = ioctl( ichan, MTIOCTOP, (char *)&magop) ;
	return ( rc );
}



int setupoas(ichan)
	int ichan;
{
#if ( defined(Machine4) | defined(Machine5)   )
	int iflag = 225, speed=13, word, sword=16384;
	char *dummy;
#ifdef TERMIOS
	tcgetattr(ichan, &pack);
	cfsetospeed(&pack, B9600);
	cfsetispeed(&pack, B9600);
	printf("flag: %o\n",pack.c_cflag);
	pack.c_cflag &= ~CRTSCTS;
	pack.c_cflag &= ~PARODD;
	pack.c_cflag &= ~PARENB;
	pack.c_cflag &= ~CRTSXOFF;
	printf("flag: %o\n",pack.c_cflag);
	tcsetattr(ichan, TCSANOW, &pack);
#else
	ioctl(ichan, TIOCGETP, &pack);
	pack.sg_flags &= ~(ECHO) ; 
	pack.sg_flags |=  (ANYP) ; 
	pack.sg_flags |=  (RAW) ; 
        pack.sg_ispeed=(char)speed;
        pack.sg_ospeed=(char)speed;
	ioctl(ichan, TIOCSETP, &pack);
	ioctl(ichan,TIOCLGET,&word);
	word=sword;
        ioctl(ichan,TIOCLSET,&word);
        ioctl(ichan,TIOCSDTR,0);
	ioctl(ichan, TIOCEXCL, dummy);
	ioctl(ichan, TIOCFLUSH, dummy);
#endif
	return(0);
#endif
}

sendoas(ichan,command,response,lr,iwt,lencom,lenresp)
int ichan, *lr, iwt, lencom, lenresp;
char *command, *response;
{
#if ( defined(Machine4) | defined(Machine5)  )
	int arg, ierrno, kntpr, i;
	unsigned int knt;
	int nonblocking=1, blocking=0, input=0, output=1, isblocking, wrote;
	extern int errno;
	extern int oasintc;
	char c, cr=13, nl=10, sp=32, pr=62;
	char *cpt, *wpt;
	extern void oasalrm();
	int secs,trys,try,kwt;
	ioctl(ichan, FIONBIO, &nonblocking);
	knt=1 ;
	while ( knt == 1)
	{	knt = read(ichan,&c,1) ; 
		if ( (knt == 1) && (oasintc == 1) )
		{	int iprt=c ; printf("discard: %d\n",iprt);
		}
	}
	oasichan=ichan;
	savehandler=signal(SIGALRM,oasalrm);
	kwt=iwt;
	if(kwt < 0) kwt = 605;
	secs = kwt/10;
	trys = kwt-secs*10;
	if(secs < 10) secs=10;
	if(noretry) secs=0;
	try = 0;
retry:
	errno = ioctl(ichan, FIONBIO, &blocking); isblocking=1;
        if(errno != 0) perror("sendoas FIONBIO");
	alarm((unsigned)secs);
	if( oasintc == 1 )
	{	fflush(stderr);
		printf("sent %d:",try);
		for ( wpt = command, i=0 ; i < lencom ; i++, wpt++) printf("%c",*wpt);
		printf("\n");
		fflush(stdout);
	}
	errno=0;
	wrote = write(ichan,command,lencom);
        if((wrote != lencom) || (errno != 0) )
	{	printf("sendoas: lencom=%d  wrote=%d\n",lencom,wrote);
		(void)perror("sendoas");
		exit(2);
	}
	errno=0;
	wrote = write(ichan,&cr,1); 
        if((wrote != 1) || (errno != 0) )
	{	printf("sendoas: wrote=%d\n",wrote);
		(void)perror("sendoas");
		exit(2);
	}
	cpt=response; cpt-- ;
	knt=0; kntpr=0;
	for ( i=0 ; i < 30 ; i++)
	{
		ierrno=0;
		errno=0;
		while(ierrno == 0)
		{	int iread;
			iread=read(ichan,&c,1);
			if( iread == 0 & isblocking == 0 ) errno=EWOULDBLOCK;
			ierrno=errno;
			if((iread == 1) )
			{	
				ierrno=errno;
				if( c != cr )
				{	cpt++ ; knt++ ; *cpt = c ;
					if(c == pr) kntpr--;
					else if(c=='<') kntpr++;
				}
			}
			else if(isblocking)
			{
				ierrno=errno;
				try++;
				if((try+1) == trys) signal(SIGALRM,savehandler);
				if( try < trys ) goto retry ;
	
			}
			else
			{
				ierrno=errno;
			}
			if( isblocking ) { ioctl(ichan, FIONBIO, &nonblocking); isblocking = 0; }
		}
		if( (ierrno == EWOULDBLOCK ) &  kntpr == -1 ) break;
		if( (ierrno != 0)  ) {  usleep((unsigned)64000) ; } 
	}
	if( knt == 0 ) knt=-1;
	else for (   ; ( (*cpt==cr) || (*cpt==sp) || (*cpt==nl) || (*cpt==pr) ) & knt > 0 ; cpt--, knt--) ;
	if( oasintc == 1 )
	{	printf("allr: knt= %d:",knt);
		for ( wpt = response, i=0 ; i < knt ; i++, wpt++) 
		{ int iprt = *wpt ; printf("AAA %c %d ",*wpt,iprt);
		}
		printf("\n");
		fflush(stdout);
	}
	*lr=knt;
	if( oasintc == 1 ) printf("returning: knt = %d\n",knt);
	if( oasintc == 1 )
	{	printf("recd:");
		for ( wpt = response, i=0 ; i < knt ; i++, wpt++) printf("%c",*wpt);
		printf("\n");
		fflush(stdout);
	}
	alarm((unsigned)0);
	signal(SIGALRM,savehandler);
#endif
}

void oasalrm(sig,code)
int sig, code;
{
#if ( defined(Machine4) | defined(Machine5)  )
	int nonblocking=1;
	printf("oasalrm: sig= %d\n",sig);
	fflush(stdout);
	ioctl(oasichan, FIONBIO, &nonblocking);
/*      If you want to re-try sendoas after a timeout omit the following line */
/*	savehandler(sig,code); */
#endif
}

statf(ichan, size, istat)
int  ichan, *size, *istat;
{
	extern int      fstat();
	int             dummy;
	dummy = fstat(ichan, &buf);
	*size =  buf.st_size;
	*istat = buf.st_mode;
}

nousleep(a)
unsigned int *a;
{
}
#endif
