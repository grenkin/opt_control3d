#ifndef INPUT_DATA_H_INCLUDED
#define INPUT_DATA_H_INCLUDED

typedef enum { COST_FUNC_J1, COST_FUNC_J2 } cost_func_t;
typedef enum { INIT_GUESS_UMIN, INIT_GUESS_UMAX } init_guess_t;

extern double ll;
extern double tt;
extern double a;
extern double alpha;
extern double kappaa;
extern double b;
extern double beta;
extern double thetab;
extern double thetainit;
extern double umin;
extern double umax;
extern double thetad;
extern cost_func_t cost_func;
extern init_guess_t init_guess_type;
extern int N;
extern int tnum;

void get_input_data();

#endif // INPUT_DATA_H_INCLUDED
