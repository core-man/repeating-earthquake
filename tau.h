# ifndef LIBTAU_H
# define LIBTAU_H

void brnset(const char *branch);
void depset(float zs);
void trtm(float delta, int *pn, float *tt, float *ray_p, float *dtdd, float *dtdh, float *dddp, char **phnm);
void emdlv(float r, float *vp, float *vs);
int tabin(const char *);
int tabout();

# endif
