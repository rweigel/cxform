/*
** cxform.h  --  prototypes and headers for Ed's coordinate transform package
*/

typedef double Vec[3];
typedef double Mat[3][3];

int cxform(const char *from,const char *to,const double et,Vec v_in,Vec v_out);

char *cxform_err(void);

/*
** Matrix multiplication and transposition
*/
void mat_transpose(Mat m_in, Mat m_out);
void mat_times_mat(Mat m1,   Mat m2, Mat m_out);
void mat_times_vec(Mat m1,   Vec v1, Vec v_out);

enum direction { FORWARD, BACK };
typedef enum direction Direction;
