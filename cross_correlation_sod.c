/******************************************************************
 * This code is used to do cross-correlation between sacfile lists
 *
 * History
 *      -     Li Sun         initial coding
 *   10/2013  Jiayuan Yao    re-organize the structures
 *   07/2015  Jiayuan Yao    set reading data error
 ******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include "complex.h"
#include <math.h>
#include "fftw3.h"
#include "sac.h"
#include "sac1.h"
#include "tau.h"
#define FL 256
#define REARTH 6371
#define eps 1.0e-5

typedef struct _ccvalue{
    int imax;
    double ccmax;
    double scale;
}ccvalue;

double htoe(int year, int mon, int day, int hour, int min, double sec);
int taper_cos(float *tpr, int n, double per);
int taper_hanning(double *tpr, int n, double per);
double lagr(double x, float a[], int n, double x0, double dx);
double dist(double evla, double evlo, double stla, double stlo);
int rsachead(int *P_flag, char **str, int len, int len1, int *len1_valid, double tk1, double tk2, double delta, char ph[]);
int rsacdata(int *P_flag,	char **str, char **str1, float **data, int *nd, double *ts, int *npts, double *b, double delta, double *tt, double *ela, double *elo, double *evdp, double *sla, double *slo, double *gcar, int *len1_valid, int len, int len1, double tk1, double tk2, char ph[], int BpFlag);
ccvalue crosscorrelation(int nd, int nd1, float *data, float *data1, double *cc);
int Getname(char str[FL], char *name);



int main(int argc, char *argv[]){

    char filen[FL], **str, **str1;
    int *npts;
    double *ela, *elo, *evdp, *sla, *slo, *gcar;
    double delta_date;

    int	*nd, *P_flag, len_tmp, i_tmp;
    char cp_tmp[FL], delims[]=".", *split_cp=NULL;

    int len = 0, i, len1, *len1_valid, len1_valid_tmp = 0, len2=0, j, len_tmp_valid = 0, num, k;
    int opt, ii, imax;

    float **data, *data_tmp;
    double *b, *tt, *ts, tk1 = 0.0, tk2 = 0.0, delta=0.025;
    float *ccf;
    char ph[8], ht[8];
    int BpFlag;

    char knetwk[20], kstnm[20], kloc[20], kcmp[20], kdate[30];
    char knetwk1[20], kstnm1[20], kloc1[20], kcmp1[20], kdate1[30];
    char fcrr_sac[FL], ccfname[FL], *cp;
    double edist;
    double inv_nfft, pw, pw1, scale, cc0, ccmax, *cc;
    struct sac_head hdr = sac_null;
    FILE *fp, *fcr;
    ccvalue ccv;


    /* read arguments */
    if (argc < 2) {
        fprintf(stderr, "Usage: %s [-t t1,t2 -d delta] masterlist lists\n", argv[0]);
        exit(1);
    }

    while ((opt = getopt(argc, argv, "t:d:p:h:b:")) != -1) {
        switch(opt) {
        case 't':
            sscanf(optarg, "%lf,%lf", &tk1, &tk2);
            break;
        case 'd':
            sscanf(optarg, "%lf", &delta);
            break;
        case 'p':
            sscanf(optarg, "%s", ph);
            break;
        case 'h':
            sscanf(optarg, "%s", ht);
            break;
        case 'b':
            sscanf(optarg, "%d", &BpFlag);
            break;
        default:
            fprintf(stderr, "Usage: %s [-t t1,t2 -d delta -p P] "
		    "list1 list2\n", argv[0]);
            exit(1);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "expect argument after the option\n");
        exit(1);
    }


    /* find line number of all lists */
    i=0;
    for (j = optind; j < argc; j++) {
        if((fp = fopen(argv[j], "r")) == NULL){
            fprintf(stderr, "open %s failed.\n", argv[j]);
            exit(1);
        }
        while(fgets(filen, sizeof(filen), fp) !=NULL){
            i++;
        }
        fclose(fp);
    }
    len=i;

    /* initialize reference phase */
    P_flag=(int*)malloc(sizeof(int)*len);
    for(i=0;i<len;i++){
        P_flag[i]=0;
    }

    /* set sacfile lists */
    str1=(char**)malloc(sizeof(char*)*len);
    for(i=0;i<len;i++){
        str1[i]=(char*)malloc(sizeof(char)*FL);
    }

    /* read main list */
    i = 0;
    if((fp = fopen(argv[optind], "r")) == NULL){
        fprintf(stderr, "open %s failed.\n", argv[optind]);
        exit(1);
    }
    while(fgets(filen, sizeof(filen), fp) !=NULL){
        filen[strlen(filen)-1] = '\0';
        strcpy(str1[i], filen);
        //printf("%s\n", str1[i]);
        i++;
    }
    fclose(fp);
    len1=i;

	/* read following list */
    for(j = optind+1; j < argc; j++){
        if((fp = fopen(argv[j], "r")) == NULL){
            fprintf(stderr, "open %s failed.\n", argv[j]);
            exit(1);
        }
        while(fgets(filen, sizeof(filen), fp) !=NULL){
            filen[strlen(filen)-1] = '\0';
            strcpy(str1[i], filen);
            //printf("%s\n", str1[i]);
            i++;
        }
        fclose(fp);
    }

    /* check valid line number */
    /* len        : number of all lines
     * len1       : line number of main list
     * len_tmp    : valid number of all lines
     * len1_valid : valid line number of main list
     */
    len1_valid=&len2;
    //fprintf(stderr,"%d	%d	%d\n",len, len1, *len1_valid);
    //fprintf(stderr,"%s:Read sachead to determine which data are chosen to read\n",argv[optind]);
    len_tmp=rsachead(P_flag,str1,len,len1,len1_valid,tk1,tk2,delta,ht);
    len1_valid_tmp = *len1_valid;
    //fprintf(stderr,"%d	%d	%d	%d\n",len, len1, *len1_valid, len_tmp);


    /* initialize variables */
    str = (char**)malloc(sizeof(char*)*len_tmp);
    for(i=0;i<len_tmp;i++){
        str[i] = (char*)malloc(sizeof(char)*FL);
    }
    data=(float**)malloc(sizeof(float*)*len_tmp);
    npts=(int*)malloc(sizeof(int)*len_tmp);
    //b=(float*)malloc(sizeof(float)*len_tmp);
    //tt=(float*)malloc(sizeof(float)*len_tmp);
    b=(double*)malloc(sizeof(double)*len_tmp);
    tt=(double*)malloc(sizeof(double)*len_tmp);
    //delta=(float*)malloc(sizeof(float)*len_tmp);
    //delta=(double*)malloc(sizeof(double)*len_tmp);

    ela=(double*)malloc(sizeof(double)*len_tmp);
    elo=(double*)malloc(sizeof(double)*len_tmp);
    sla=(double*)malloc(sizeof(double)*len_tmp);
    slo=(double*)malloc(sizeof(double)*len_tmp);
    gcar=(double*)malloc(sizeof(double)*len_tmp);
    evdp=(double*)malloc(sizeof(double)*len_tmp);

    //ts=(float*)malloc(sizeof(float)*len_tmp);
    ts=(double*)malloc(sizeof(double)*len_tmp);
    nd=(int*)malloc(sizeof(int)*len_tmp);

    /* read sac data */
    //fprintf(stderr,"%s:Read sacdata\n",argv[optind]);
    len_tmp_valid = rsacdata(P_flag,str,str1,data,nd,ts,npts,b,delta,tt,ela,elo,evdp,sla,slo,gcar,len1_valid,len,len1,tk1,tk2,ht,BpFlag);

    free(P_flag);
    for (i=0; i < len; i++) {
        free(str1[i]);
    }
    free(str1);

    /* check valid sacfile number */
    //fprintf(stderr,"len: %d		len_valid: %d\n",len_tmp, len_tmp_valid);
    //fprintf(stderr,"len1_all: %d		len1_valid: %d\n",len1_valid_tmp, *len1_valid);
    if((len_tmp     != len_tmp_valid) ||
       (*len1_valid != len1_valid_tmp)){
        fprintf(stderr,"**********The valid datalist is not the same.***********\n");
        fprintf(stderr,"******************%d	%d.**************\n", len_tmp, len_tmp_valid);
        exit(1);
    }

    /* output file name */
    snprintf(ccfname, 256, "%s.%s.ccor", argv[optind],ph);
    if ((fcr = fopen(ccfname, "w")) == NULL) {
        fprintf(stderr, "Failed to open ccor\n");
        exit(1);
    }

    /****************************************************
    ****    do cross-correlate in frequency domain   ****
    ****************************************************/
    //fprintf(stderr,"cross-correlate\n");
    //for(i = 0; i < len_tmp; i++){
    for(i = 0; i < *len1_valid; i++){
        //fprintf(stderr,"%s\n", str[i]);

        cp = strrchr(str[i], '/');
        if (cp == NULL) {
            cp = str[i];
        }
        else {
            cp++;
        }

        /* get network, station, location etc. */
        /* you may change them for you directory structure and file name */
        //sscanf(cp, "%[^.].%[^.].%[^.].%[^.].%[^.]", knetwk, kstnm, kloc,kcmp, kdate);
        //sscanf(cp, "%[^.].%[^.].%[^.].%[^.].", kdate, knetwk, kstnm,kloc);
        sscanf(cp, "%[^.].%[^.].%[^.].", knetwk, kstnm, kloc);
        Getname(str[i], kdate);
        //fprintf(stderr,"%s %s %s %s\n", kdate, knetwk, kstnm, kloc);

        strcpy(cp_tmp,cp);
        split_cp=strtok(cp_tmp,delims);
        i_tmp=0;
        while(split_cp!=NULL){
            //if(strcmp(split_cp,"BHZ")==0 && i_tmp==4){
            if(strcmp(split_cp,"BHZ") == 0 && i_tmp==2){
                strcpy(kloc,"");
                //fprintf(stderr,"############%s\n",split_cp);
            }
            i_tmp++;
            split_cp=strtok(NULL,delims);
        }

        //fprintf(stderr, "cross correlation %d	%d %s\n",i,len_tmp, kdate);

        //for(j = i+1;j < len_tmp; j++){
        for(j = i+1;j < len_tmp_valid; j++){
            //edist = dist(ela[i], elo[i], ela[j], elo[j]);
            //fprintf(stderr,"edist is %f\n",edist);
            edist = (REARTH-evdp[i])*(REARTH-evdp[i])+(REARTH-evdp[j])*(REARTH-evdp[j])-2*(REARTH-evdp[i])*(REARTH-evdp[j])*cos(dist(ela[i], elo[i], ela[j], elo[j]));
            edist= sqrt(edist);
            //fprintf(stderr,"%s\n", str[j]);
            //fprintf(stderr,"edist is %f\n",edist);
            //fprintf(stderr,"edist is %f	%f\n",evdp[i], evdp[j]);

            cp = strrchr(str[j], '/');
            if (cp == NULL) {
                cp = str[j];
            }
            else {
                cp++;
            }

            //sscanf(cp, "%[^.].%[^.].%[^.].%[^.].%[^.]", knetwk1, kstnm1, kloc1,kcmp1, kdate1);
            sscanf(cp, "%[^.].%[^.].%[^.].", knetwk1, kstnm1, kloc1);
            Getname(str[j], kdate1);
            //sscanf(cp, "%[^.].%[^.].%[^.].%[^.].", kdate1, knetwk1, kstnm1,kloc1);

            strcpy(cp_tmp,cp);
            //fprintf(stderr,"%s\n",cp_tmp);
            split_cp=strtok(cp_tmp,delims);
            i_tmp=0;
            while(split_cp!=NULL){
                //if(strcmp(split_cp,"BHZ")==0 && i_tmp==4){
                if(strcmp(split_cp,"BHZ")==0 && i_tmp==2){
                    strcpy(kloc1,"");
                    //fprintf(stderr,"hahahahah%s\n",split_cp);
                }
                i_tmp++;
                split_cp=strtok(NULL,delims);
            }

            delta_date=fabs(atol(kdate)-atol(kdate1));
            //fprintf(stderr, "cross correlation %ld	%ld\n",atol(kdate),atol(kdate1));
            //fprintf(stderr, "cross correlation %d	%d	%s %s %f	%f\n",i, j, kdate,kdate1,edist,delta_date);




            //if(edist <= 60 && delta_date >= 3000 && (strcmp(knetwk,knetwk1)==0) && (strcmp(kstnm,kstnm1)==0)){ //The two earthquakes must be close between hypocenter and far from origin time.
            //if(edist <= 60 && delta_date > 0 && (strcmp(knetwk,knetwk1)==0) && (strcmp(kstnm,kstnm1)==0)){

            /*relocation*/
            //if(edist <= 60 && delta_date > 0){

            /*calculate ccf*/
            if(1){
                //fprintf(stderr, "cross correlation %s %s	%f\n",kdate,kdate1,edist);
                if ((cc = malloc((nd[i] + nd[j] - 1) * sizeof(*cc))) == NULL) {
                    fprintf(stderr, "allocation failed for data, npts\n");
                    continue;
                }
                if((ccf = malloc((nd[i] + nd[j] - 1) * sizeof(*ccf))) == NULL) {
                    fprintf(stderr, "allocation failed for data, npts\n");
                    continue;
                }

                ccv=crosscorrelation(nd[i], nd[j], data[i], data[j], cc);

                //ccmax=ccv.ccmax; imax=ccv.imax; scale=ccv.scale;
                //fprintf(stderr,"%20.18f\n",delta;

                if ((fabs(ccv.ccmax*ccv.scale)) > 0.0) {
					fprintf(fcr,"%s  %s  %f  %f\n",str[i], str[j],ccv.ccmax*ccv.scale, -(ccv.imax*delta+ts[i]-ts[j]));
                    //fprintf(fcr,"%s  %s  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",str[i],str[j],ccv.ccmax*ccv.scale,-(ccv.imax*delta+b[i]-b[j]),edist,sla[i],slo[i],gcar[i],gcar[j],ela[i],elo[i],ela[j],elo[j]);
                    //fprintf(fcr,"%s  %s  %f  %f  %f\n",str[i], str[j],ccv.ccmax*ccv.scale,-(ccv.imax*delta+ts[i]-ts[j]),edist);
                    //fprintf(fcr,"%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",ccv.ccmax*ccv.scale,-(ccv.imax*delta+ts[i]-ts[j]),edist,sla[i],slo[i],gcar[i],gcar[j],ela[i],elo[i],ela[j],elo[j]);
                    //fprintf(fcr,"%s  %s  %f  %f  %f %f %f\n",str[i], str[j],ccv.ccmax*ccv.scale,-ccv.imax*delta,-(ts[i]-ts[j]),-(ccv.imax*delta+ts[i]-ts[j]),edist);
                    //fprintf(fcr,"%s  %s  %f  %f  %f %f\n",str[i],str[j],ccv.ccmax*ccv.scale,-ccv.imax*delta,-(ts[i]-ts[j]),edist);
                }

                //if (((fabs(ccv.ccmax * ccv.scale)) > 1.0) && ((fabs(ccv.imax * delta)) < 20)) {
                /*if ((fabs(ccv.ccmax * ccv.scale)) >= 0.8) {
                    //snprintf(fcrr_sac, 256,"/home/yaojy/Test/%s.%s.%s.%s.%s.%s.cr.SAC", knetwk, kstnm,kloc,kloc1, kdate, kdate1);
                    //snprintf(fcrr_sac, 256,"/public/home/yaojy/doublet_search/CC/%s.%s.%s.%s.%s.%s.cr.SAC", knetwk, kstnm,kloc,kloc1, kdate, kdate1);
                    if ((fp = fopen(fcrr_sac, "w")) == NULL) {
                        fprintf(stderr, "Failed to open %s\n", fcrr_sac);
                        exit(1);
                    }

                    for(k = 0; k < (nd[i] + nd[j]) - 1; k++)
                        ccf[k] = cc[k] * ccv.scale;

                    hdr.t1 = -12345.;
                    hdr.nzyear = -12345;
                    hdr.nzjday = -12345;
                    hdr.nzhour = -12345;
                    hdr.nzmin = -12345;
                    hdr.nzsec = -12345;
                    hdr.nzmsec = -12345;
                    hdr.stla = sla[i];
                    hdr.stlo = slo[i];
                    strncpy(hdr.kstnm, kstnm, 8);
                    strncpy(hdr.knetwk, knetwk, 8);

                    if ((tk1 + tk2) > 0)
                        hdr.b = (-nd[i] + 1) * delta + (ts[j] - ts[i]);
                    else
                        hdr.b = (-nd[i] + 1) * delta + (b[j] - b[i]);
                    hdr.npts = nd[i] + nd[j] - 1;
                    hdr.e = hdr.b + (hdr.npts - 1) * delta;
                    hdr.o = 0.0;
                    fwrite(&hdr, sizeof(hdr), 1, fp);
                    fwrite(ccf, sizeof(*ccf), nd[i] + nd[j] - 1, fp);
                    fclose(fp);
                }*/

                free(cc);
                free(ccf);
            } /* if */

        } /* for j */

        free(str[i]);
        free(data[i]);
    } /* for i */

    /* free memory */
    for(i = *len1_valid;i < len_tmp; i++){
        free(str[i]);free(data[i]);
    }
    free(str);free(data);
    free(ts);free(nd);
    free(npts);free(b);free(tt);
    //free(delta);
    free(sla);free(slo);free(ela);free(elo);free(gcar);free(evdp);
    fclose(fcr);

}



/*************
 * functions
 *************/
double htoe(int year, int mon, int day, int hour, int min, double sec)
{
    double e;
    int i, days = 0;
    int d[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    if (year > 1970) {
        for (i = 1970; i < year; i++) {
            days += 365;
            if ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0)
                days++;
        }
    } else if (year < 1970) {
        for (i = year; i < 1970; i++) {
            days -= 365;
            if ((i % 4 == 0 && i % 100 != 0) || i % 400 == 0)
                days--;
        }
    }

    if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0)
        d[1]++;

    if (mon < 0)
        days += day - 1;
    else {
        for (i = 0; i < mon - 1; i++)
            days += d[i];
        days += day - 1;
    }

    e = days * 86400.0;
    e += hour * 3600.0 + min * 60.0 + sec;

    return e;
}


int taper_cos(float *tpr, int n, double per)
{
    int i, nptr;
    double omega, v;

    if (per > 0.5)
        return -1;

    for (i = 0; i < n; i++)
        tpr[i] = 1.0;

    nptr = (int) ((n + 1) * per);
    if (nptr < 2)
        nptr = 2;
    omega = 0.5 * M_PI / nptr;

    for (i = 0; i < nptr; i++) {
        v = sin(omega * i);
        tpr[i] *= v;
        tpr[n - 1 - i] *= v;
    }

    return 0;
}


int taper_hanning(double *tpr, int n, double per)
{
    int i, nptr;
    double omega, v, f0 = 0.5, f1 = 0.5;

    if (per > 0.5)
        return -1;

    for (i = 0; i < n; i++)
        tpr[i] = 1.0;

    nptr = (int) ((n + 1) * per);
    if (nptr < 2)
        nptr = 2;
    omega = M_PI / nptr;

    for (i = 0; i < nptr; i++) {
        v = f0 - f1 * cos(omega * i);
        tpr[i] *= v;
        tpr[n - 1 - i] *= v;
    }

    return 0;
}


double lagr(double x, float a[], int n, double x0, double dx)
{
    int n0;
    double d12, d01, dm0, dm1, d02, dm2, xn0, t, t1, t2;
    double y;

    n0 = (int)((x - x0) / dx);
    if (n0 < 1) n0 = 1;
    if (n0 > n - 3) n0 = n - 3;
    d12 = a[n0 + 2] - a[n0 + 1];
    d01 = a[n0 + 1] - a[n0];
    dm0 = a[n0] - a[n0 - 1];

    dm1 = d01 - dm0;
    d02 = d12 - d01;
    dm2 = d02 - dm1;

    xn0 = x0 + n0 * dx;
    t = (x - xn0) / dx;
    t1 = t * (t - 1.0) / 2.0;
    t2 = t1 * (t + 1.0) / 3.0;

    y = a[n0] + t * d01 + t1 * dm1 + t2 * dm2;

    return y;
}

double dist(double evla, double evlo, double stla, double stlo)
{
    double delta;
    double s1, s2, c1, c2, cd, lon_dist;

    evla *= M_PI / 180.0;
    stla *= M_PI / 180.0;
    lon_dist = (evlo - stlo) * M_PI / 180.0;

    /* Calculate geocentric latitude.
     * Use GRS80 reference ellipsoid.
     * a = 6378137.0 meters (equatorial radius)
     * b = 6356752.0 meters (polar radius)
     * f = (b/a)^2 = 0.99330552180201
     */
    s1 = 0.99330552180201 * sin(evla);
    s2 = 0.99330552180201 * sin(stla);
    evla = atan2(cos(evla), s1);
    stla = atan2(cos(stla), s2);

    s1 = sin(evla);
    s2 = sin(stla);
    c1 = cos(evla);
    c2 = cos(stla);
    cd = cos(lon_dist);

//    delta = acos(c1 * c2 + s1 * s2 * cd) / M_PI * 180.0;
    delta = acos(c1 * c2 + s1 * s2 * cd);

    return delta;
}


/* check valid line number */
int rsachead(int *P_flag, char **str, int len, int len1, int *len1_valid, double tk1, double tk2, double delta, char ph[]){
    int i, num, nd, ns, npts, len_valid;
    float b, ts, tt;
    struct sac_head hdr = sac_null;
    FILE *fp;

    len_valid=0;
    for(i=0;i<len;i++){
        //fprintf(stderr,"%s\n",str[i]);
        if((fp = fopen(str[i], "r")) == NULL){
            fprintf(stderr, "## open %s failed.\n", str[i]);
            exit(1);
            continue;
        }
        if(fread(&hdr, sizeof(hdr), 1, fp) != 1){
            fprintf(stderr, "read failed for header of %s\n",str[i]);
            continue;
        }
        fclose(fp);

        npts = hdr.npts;
        b = hdr.b;
        //delta = hdr.delta;
        //tt = hdr.t1;
        //fprintf(stderr,"%s\n",ph);

        /* set reference time */
        if(strcmp(ph,"t1")==0){
            tt = hdr.t1;
            //fprintf(stderr,"****%s	%f\n", ph, tt);
        }else if(strcmp(ph,"t2")==0){
            tt = hdr.t2;
        }else if(strcmp(ph,"t3")==0){
            tt = hdr.t3;
        }else if(strcmp(ph,"t4")==0){
            tt = hdr.t4;
        }else if(strcmp(ph,"t5")==0){
            tt = hdr.t5;
        }else if(strcmp(ph,"t6")==0){
            tt = hdr.t6;
        }else if(strcmp(ph,"t7")==0){
            //tt = hdr.t6;
            tt = hdr.t7;
        }else if(strcmp(ph,"t8")==0){
            tt = hdr.t8;
        }else if(strcmp(ph,"t9")==0){
            tt = hdr.t9;
        }else if(strcmp(ph,"t0")==0){
            tt = hdr.t0;
        }else {
            tt = hdr.o;
            //fprintf(stderr,"%f\n", tt);
        }
        //tt = hdr.t2;
        //fprintf(stderr,"****%s	%f****\n", ph, tt);

        if(fabs(tt+12345)< 1.e-5){continue;}

        /* find start point and point number */
        ns=0;
        if ( (tk2+tk1)> 0 ) {
            nd=(tk1+tk2+0.00001)/delta;

            for(num = 0; num < npts; num++){
                ts = b + num*delta;
                if (ts < (tt - tk1))
                    ns++;
                else
                    break;
            }

            if (num==0) {
                nd=(tt+tk2-ts+0.00001)/delta+1;
            }
            //fprintf(stderr,"%d  %d  %d\n",nd,npts,ns);

            if (npts < (ns + nd)) {
                nd=npts-ns;
            }
        } else {
            nd = npts;
            ts=b;
        }

        //fprintf(stderr,"%s\n",str[i]);
        //fprintf(stderr,"The master event:%20.15f %20.18f  %20.18f  %d\n",ts, b,delta,num);

        if(nd < 2 || ts >= tt + tk2){
            //fprintf(stderr,"%d	%f	%f\n%s\n", nd, ts, tt, str[i]);
            continue;
        } // The time window is not in the time interval of data.

        //fprintf(stderr,"%d	%f	%f\n%s\n", nd, ts, tt, str[i]);
        len_valid++;

        if (i < len1) {
            *len1_valid=len_valid;
        }
        P_flag[i]=1;
    }

    //fprintf(stderr,"%d %d\n",len_valid, *len1_valid);
    return len_valid;
}


/* read sac data */
int rsacdata(int *P_flag, char **str, char **str1, float **data, int *nd, double *ts, int *npts, double *b, double delta, double *tt, double *ela, double *elo, double *evdp, double *sla, double *slo, double *gcar, int *len1_valid, int len, int len1, double tk1, double tk2, char ph[], int BpFlag){

    int	i, j, k, num, ns;
    float *data_tmp;
    struct sac_head hdr = sac_null;
    FILE *fp, *fp1;

    // for bandpass
    double low=0.8, high=2;
    //double low=0.5, high=1;
    //double attenuation=0, transition_bandwidth=0;
    double attenuation=30, transition_bandwidth=0.3;
    //int nlen;
    int order=2,passes=1;
    //xapiir(yarray,nlen,SAC_BUTTERWORTH,transition_bandwidth,attenuation,order,SAC_BANDPASS,low,high,delta_d,passes);


    /*if((fp1 = fopen("filelist", "w")) == NULL){
        fprintf(stderr, "****************open failed.\n");
        exit(1);
    }*/


    k=0;
    for(i=0;i<len;i++){
        //fprintf(stderr,"%d	%d\n",P_flag[i], len);
        if(P_flag[i] == 1){
            if((fp = fopen(str1[i], "r")) == NULL){
                fprintf(stderr, "****************open %s failed.\n", str1[i]);
                continue;
            }
            if(fread(&hdr, sizeof(hdr), 1, fp) != 1){
                fprintf(stderr, "read failed for header of %s\n",str1[i]);
                continue;
            }



			strcpy(str[k], str1[i]);
//			fprintf(stderr,"***********###############%s	%s\n",str[k],str1[i]);
//			fprintf(fp1,"%s\n",str[k]);
			npts[k] = hdr.npts;
			b[k] = hdr.b;
//			delta[k] = hdr.delta;
//			tt[k] = hdr.t1;
			ela[k] = hdr.evla;
			elo[k] = hdr.evlo;
			sla[k] = hdr.stla;
			slo[k] = hdr.stlo;
			gcar[k] = hdr.gcarc;
			//evdp[k] = hdr.evdp;
			evdp[k] = hdr.evdp / 1000.0;	// the unit is meter for sod
//		fprintf(stderr,"%s\n",ph);

		if(strcmp(ph,"t1")==0){
			tt[k] = hdr.t1;
		}else if(strcmp(ph,"t2")==0){
			tt[k] = hdr.t2;
		}else if(strcmp(ph,"t3")==0){
			tt[k] = hdr.t3;
		}else if(strcmp(ph,"t4")==0){
			tt[k] = hdr.t4;
		}else if(strcmp(ph,"t5")==0){
			tt[k] = hdr.t5;
		}else if(strcmp(ph,"t6")==0){
			tt[k] = hdr.t6;
		}else if(strcmp(ph,"t7")==0){
//			tt[k] = hdr.t6;
			tt[k] = hdr.t7;
		}else if(strcmp(ph,"t8")==0){
			tt[k] = hdr.t8;
		}else if(strcmp(ph,"t9")==0){
			tt[k] = hdr.t9;
		}else if(strcmp(ph,"t0")==0){
			tt[k] = hdr.t0;
		}else {
			tt[k] = hdr.o;
		}

			if(fabs(tt[k]+12345)< 1.e-5){continue;}
			ns=0;

			if((tk2+tk1)>0) {
    			nd[k]=(tk1+tk2+0.00001)/delta;
    			for(num = 0; num < npts[k]; num++){
    				ts[k]=b[k]+num*delta;
      				if(ts[k] < (tt[k] - tk1))
						ns++;
      				else
      					break;
    			}
    			if(num == 0) {
    	  			nd[k]=(tt[k]+tk2-ts[k]+0.00001)/delta+1;
    			}
	  			if(npts[k] < (ns + nd[k]))
    				nd[k]=npts[k]-ns;
			} else {
    			nd[k] = npts[k];
				ts[k]=b[k];
			}
			if(nd[k] < 2 || ts[k] >= tt[k] + tk2){continue;}

			if((data_tmp = (float*)malloc(npts[k] * sizeof(float))) == NULL){
    			fprintf(stderr, "allocation failed for data %d.\n",i);
      			continue;
    		}
			if(fread(data_tmp, sizeof(float), npts[k], fp) != npts[k]){
    			fprintf(stderr, "read failed for data of %s\n",str1[i]);
    	  		continue;
    		}

			if ( BpFlag == 1 ) {	// do the bandpass
	    		xapiir(data_tmp, npts[k], SAC_BUTTERWORTH,transition_bandwidth, attenuation,order,SAC_BANDPASS,low, high,delta, passes);
			}

			if((data[k] = (float*)malloc(nd[k]*sizeof(float))) == NULL){
    			fprintf(stderr, "Allocation failed for data %d.\n",i);
    	  		continue;
    		}
			for(j=0;j<nd[k];j++){
				data[k][j]=data_tmp[ns+j];
			}


			k++;
			if(i < len1){*len1_valid = k;}


			free(data_tmp);
			fclose(fp);

		}
//	  	fprintf(stderr,"The master event:%d	%d\n",i,k);
//		free(str1[i]);
	}
//	fclose(fp1);
//	fprintf(stderr,"%d	%d	%d\n", k,i-1, len);
	return k;

}


ccvalue crosscorrelation(int nd, int nd1, float *data, float *data1, double *cc){
	int k, num, nw, nfft;
	int imax;
	double scale, ccmax;
  double *d, *d1, *dtap;
  double inv_nfft, pw, pw1, cc0;
  double complex *fft_in, *fft_out, *fft_out1;
  fftw_plan pf, pb;
	ccvalue ccv;

	if(nd >= nd1) {
		for(k = 2 * nd, nfft = 1; nfft < k; nfft += nfft);
	  	nw = nd;
	}
	else{
		for(k = 2 * nd1, nfft = 1; nfft < k; nfft += nfft);
	  	nw = nd1;
	}
	inv_nfft = 1.0 / (double) nfft;

	if((d = malloc(3 * nw * sizeof(*d))) == NULL) {
  	fprintf(stderr, "allocation failed for d.\n");
    exit(1);
  }
	d1 = d + nw;
  dtap = d1 + nw;
  taper_hanning(dtap, nw, 0.05);

//for (k = ns, num = 0; k < npts[i], num < nd; k++){
  for ( num = 0; num < nd; num++){
//	d[num] = data[k] * dtap[num];
		d[num] = data[num] * dtap[num];
//  num++;
  }

	for (pw = 0.0, num = 0; num < nd; num++)
	{
  	pw += d[num] * d[num];
//		fprintf(stderr,"%f	%f\n",pw,d[num]*d[num]);
  }
//for (k = ns1, num = 0; k < npts1[j], num < nd1; k++) {
	for (num = 0; num < nd1; num++) {
//	d1[num] = data1[k] * dtap[num];
		d1[num] = data1[num] * dtap[num];
//	num++;
	}

	for(pw1 = 0.0, num = 0; num < nd1; num++)
		pw1 += d1[num] * d1[num];

	scale = 1.0 / sqrt(pw * pw1);
//	fprintf(stderr,"%d %f %f %20.18f\n",nw,pw,pw1,scale);

	if((fft_in = malloc(3 * nfft * sizeof(*fft_in))) == NULL){
		fprintf(stderr, "allocation failed for fft_in\n");
		exit(1);
	}

	fft_out = fft_in + nfft;
	fft_out1 = fft_out + nfft;
	pf = fftw_plan_dft_1d(nfft, fft_in, fft_out, FFTW_FORWARD,FFTW_ESTIMATE);
	pb = fftw_plan_dft_1d(nfft, fft_in, fft_out1, FFTW_BACKWARD,FFTW_ESTIMATE);

	for (k = 0; k < nfft; k++)
		if(k < nd)
			fft_in[k] = d[k];
		else
			fft_in[k] = 0.0;
	fftw_execute_dft(pf, fft_in, fft_out);
  for (k = 0; k < nfft; k++)
    if (k < nd1)
    	fft_in[k] = d1[k];
    else
      fft_in[k] = 0.0;
	fftw_execute_dft(pf, fft_in, fft_out1);
	for (k = 0; k < nfft; k++)
		fft_out[k] *= conj(fft_out1[k]) * inv_nfft;
	fftw_execute_dft(pb, fft_out, fft_in);
  fftw_destroy_plan(pf);
  fftw_destroy_plan(pb);

	imax = 0;
	ccmax = creal(fft_in[0]);
	cc0 = fabs(ccmax);
	for (num = 0, k = nd-1; num < nd; num++){
		cc[k] = creal(fft_in[num]);
	  if (fabs(cc[k]) > cc0){
			imax = num;
	  	cc0 = fabs(cc[k]);
	  	ccmax = cc[k];
	  }
	  k--;
	}

	for (num = 1, k = nd; num < nd1 ; num++) {
		cc[k] = creal(fft_in[nfft - num]);
		if (fabs(cc[k]) > cc0) {
			imax = -num;
			cc0 = fabs(cc[k]);
		  ccmax = cc[k];
	  }
	  k++;
  }
	free(d);free(fft_in);
	ccv.imax=imax;
	ccv.scale=scale;
	ccv.ccmax=ccmax;
	return ccv;
}


int Getname(char str[FL], char *name)
{

	char cp[FL], *split_cp = NULL, delims[] = "/";
	int i = 0;

	strcpy(cp, str);
	split_cp = strtok(cp, delims);
	//fprintf(stderr, "#####****%s*****\n", split_cp);
	//fprintf(stderr, "#####****%s*****\n", cp);

	while ( (split_cp = strtok(NULL, delims)) && i < 0) {
		strcpy(cp, split_cp);
		//fprintf(stderr, "((((****%s*****\n", cp);
		//fprintf(stderr, "((((****%s*****\n", split_cp);
		i++;
	}
	//fprintf(stderr, "((((****%s*****\n", split_cp);
	strcpy(name, split_cp);

	return 1;

}
