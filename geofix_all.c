//#define ALG_TEST
//#define VERBOSE
//#define VERBOSE1
//#define PARAM_TEST

#ifndef ALG_TEST
#include "Python.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

// 全局变量
#define TABLE_SIZE (660 * 450 * 4)

#define SCALE_FACTOR  100000.0
// (1.0 / SCALE_FACTOR)
#define INVERSE_SCALE 0.00001
// 1 / 10.0;
const double COEF_REV_10 = 0.1; 

#define NEIBOR_NUM 7
#define PARAM_NUM 3

typedef long long int64;

typedef struct _Point {
	double lo, la;
} Point;

// 保存Sogou坐标
typedef struct {
	int lo, la;
} PointInt2D;

// 一个经度对应110km, 1e5之后一个单位对应1.1m
// 采用int表示中间数据 1e5坐标系，进度为1.1m
typedef struct _Point2 {
	int lo0, la0; 
	int lo1, la1;
	struct  _Point2 *next;
} Point2;

typedef struct _Node {
	int64 dist;
	const Point2* pt;
} Node;

#define BUCKET_SIZE_LO  66
#define BUCKET_SIZE_LA  50

PointInt2D gps2Sogo[660 * 2][450 * 2];

Point2 points[TABLE_SIZE];
Point2* buckets[BUCKET_SIZE_LO][BUCKET_SIZE_LA];

// 只针对当前给定的数据
// Lo: 72.000000 - 137.900000
// La: 9.997800 - 59.400000
// 以每一个经纬度为一个单位
//const double kMaxY = 7007802;
//const double kMinY = 2039740;
//const double kMaxX = 15038536;
//const double kMinX = 8327213;

#define MIN_LO 8327213
#define MAX_LO (15038536 + 100)

#define MIN_LA  2039740
#define MAX_LA  (7007802 + 100)

static int max(int a, int b) {
    return a > b ? a : b;
}
static int min(int a, int b) {
    return a < b ? a : b;
}

#define GET_LO_IDX(lo) 		((int) ( (lo - MIN_LO) * INVERSE_SCALE))
#define GET_LA_IDX(la) 		((int) ( (la - MIN_LA) * INVERSE_SCALE))

// 函数声明
static void legs (double a[PARAM_NUM][PARAM_NUM], int n, double b1[], double x1[], 
                  double b2[], double x2[], int indx[]);

static void elgs (double a[PARAM_NUM][PARAM_NUM], int n, int indx[]);
static void initializeGoogle2Sogou(const char* path);
static Point interpolation(double lo, double la,
                           Node heap[NEIBOR_NUM]);

// Implementation
////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void initializeGps2Sogou(const char* path) {
    int lo, la;
    int col, row;
    FILE* fid;
    
    // printf("Path: %s\n", path);

	fid = fopen(path, "rt");

	int counter = 1;
	Point2 lastPts[2];
	for (col = 0; col < 450 * 2; col++) {
	    for (row = 0; row < 660 * 2; row++) {
	        int i = col * 660 * 2 + row;

            int cnt = fscanf(fid, "%d %d", &lo, &la);
            if (cnt < 2) {
                printf("Error reading data.....\n");
                exit(0);
            }
            if (i > 1) {
                lo = 2 * lastPts[counter & 0x01].lo0 - lastPts[(counter  - 1)& 0x01].lo0 - lo;
                la = 2 * lastPts[counter & 0x01].la0 - lastPts[(counter  - 1)& 0x01].la0 - la;
                counter++;
                lastPts[counter & 0x01].lo0 = lo;
                lastPts[counter & 0x01].la0 = la;
            } else {
                lastPts[i].lo0 = lo;
                lastPts[i].la0 = la;
            }
    
            PointInt2D* pt2 = &gps2Sogo[row][col];
            pt2->lo = lo;
            pt2->la = la;	        
	    }
	}

	fclose(fid);
}
static void initializeGoogle2Sogou(const char* path) {
    int  i, indexLo, indexLa;
    int lo, la;
    FILE* fid;
    
    printf("Path: %s\n", path);

	fid = fopen(path, "rt");
	memset(buckets, 0, sizeof(buckets));

	int counter = 1;
	Point2 lastPts[2];
	for (i = 0; i < TABLE_SIZE; i++) {
		int cnt = fscanf(fid, "%d %d", &lo, &la);
		if (cnt < 2) {
		    printf("Error reading data.....\n");
		    exit(0);
		}
		if (i > 1) {
		    lo = 2 * lastPts[counter & 0x01].lo0 - lastPts[(counter  - 1)& 0x01].lo0 - lo;
		    la = 2 * lastPts[counter & 0x01].la0 - lastPts[(counter  - 1)& 0x01].la0 - la;
		    counter++;
		    lastPts[counter & 0x01].lo0 = lo;
		    lastPts[counter & 0x01].la0 = la;
		} else {
		    lastPts[i].lo0 = lo;
		    lastPts[i].la0 = la;
		}

		// printf("%d %d\n", lo, la);

		Point2 * pt2 = points + i;
		
		pt2->lo0 = lo; // * INVERSE_SCALE;
		pt2->la0 = la; // * INVERSE_SCALE;

		// 必须是整数除法
		// 差值时取值为: 10000 1e5单位，及11km
		// 在[11km * 11km]的范围内进行线性差值, 如何控制误差呢?
		pt2->lo1 = (i % 1320) * 5000 + 7200000; // (i % 660) * COEF_REV_10 + 72.0;
		pt2->la1 = (i / 1320) * 5000 + 1000000; // (i / 660) * COEF_REV_10 + 10.0;
		
		//printf("PointXX: %.6f, %.6f\n", pt2->lo1 * INVERSE_SCALE, pt2->la1 * INVERSE_SCALE);

		// 按照(lo0, la0)来分桶(二纬分桶)
		indexLo = GET_LO_IDX(pt2->lo0);
		indexLa = GET_LA_IDX(pt2->la0);
		indexLa = min(max(indexLa, 0), BUCKET_SIZE_LA - 1);
	    indexLo = min(max(indexLo, 0), BUCKET_SIZE_LO - 1);

		pt2->next = buckets[indexLo][indexLa];
		buckets[indexLo][indexLa] = pt2;
	}


	fclose(fid);
	// printf("Data Loading Complete...\n");	
}

// 最后一个元素保存这新放入的元素
static void shift_up(Node* nodes, int heapLen) {
	// shift up
	int hole = heapLen - 1;
	const Node value = nodes[hole];

	while(hole > 0) {
	    int parent = (hole - 1) >> 1;
	    if (nodes[parent].dist >= value.dist) {
	        break;

	    } else {
	    	nodes[hole] = nodes[parent];
	    	hole = parent;
	    }
	}
	nodes[hole] = value;
}

// 第0个元素被删除，需要将最后一个元素从0-th的位置调整到合适的位置
// 需要维持一个最大优先级队列
static void shift_down(Node* nodes, int heapLen) {
	int hole, maxIdx, left, right;
	int64 maxDist;
	Node value;

	if (heapLen == 0) {
        return;
	}

    value = nodes[heapLen];
	hole = 0;

	left = hole * 2 + 1;
	while(left < heapLen) {
	    maxIdx = hole;
	    maxDist = value.dist;

	    if (nodes[left].dist > value.dist) {
	        maxIdx = left;
	        maxDist = nodes[left].dist;
	    }

	    // Right
	    right = left + 1;
	    if (right < heapLen) {
	        if (nodes[right].dist > maxDist) {
	        	maxIdx = right;
	        	maxDist = nodes[right].dist;
	        }
	    }

	    if (maxIdx != hole) {
	        nodes[hole] = nodes[maxIdx];
	        hole = maxIdx;
	        // For Next Iteration
	        left = hole * 2 + 1;

	    } else {
	        break;
	    }
	}
	
	nodes[hole] = value;
}

//static void printHeap(const Node*heap, int len) {
//    int i;
//    printf("Heap: ");
//    for(i = 0; i < len; i++) {
//        printf("%.6f ", heap[i].dist);
//    }
//    printf("\n");
//}
static void processList(const Point2* list, const int lo, const int la, 
                        Node* heap, int* len) {
	
	int64 errLo, errLa, dist;
	int count;
	
    if (!list) {
        return;
    }
    count = 0;
    while(list) {
    	errLo = list->lo0 - lo;
    	errLa = list->la0 - la;
    	

    	// 采用 dist^2来比较大小
    	dist = (long long)errLo * errLo + (long long)errLa * errLa;
#ifdef VERBOSE1
    	if (list->lo0 == 12957035 && list->la0 == 4838017) {
    	    printf("== %15lld, %10lld %10lld %d\n", dist, errLo, errLa, (int)sqrt(dist));
    	}
#endif
        
        if (*len < NEIBOR_NUM) {
            Node* node = heap + ((*len)++);
            node->dist = dist;
            node->pt = list;
            
            shift_up(heap, *len);
        } else {
        	if (dist < heap[0].dist) {
        	    
        	    (*len)--;
        	    shift_down(heap, *len);
        	    
        	    Node* node = heap + ((*len)++);
        	    node->dist = dist;
        	    node->pt = list;
        	    shift_up(heap, *len);
        	}
        }
        count++;
        list = list->next;
    }
    
    // printf("List Length: %d\n", count);
}
static Point _sogouToGoogleMars(const double lo, const double la) {
    // clock_t start = clock();
#ifdef VERBOSE
    printf("Input: %d, %d\n", (int)lo, (int)la);
#endif
	// 初始化数据
    Node heap[NEIBOR_NUM];
    //double buff0[N], buff1[N];
	
	const int lo1e5 = (int)lo;
	const int la1e5 = (int)la;

    int loIndex = GET_LO_IDX(lo1e5);
    int laIndex = GET_LA_IDX(la1e5);
	loIndex = min(max(loIndex, 0), BUCKET_SIZE_LO - 1);
	laIndex = min(max(laIndex, 0), BUCKET_SIZE_LA - 1);
	
	// printf( "Lo Index: %d, La Index: %d\n", loIndex, laIndex); 
	int len = 0;
	int i, j;
	int64 step = 0;


	

	// step * scale (scale == 1)
	while (1) {
	    int step1 = step - 1;
	    int to_continue = (len < NEIBOR_NUM) || (heap[0].dist > (step1 * step1 * 10000000000ll));
#ifdef VERBOSE1
        int64 step2 = (step * step * 10000000000ll);
        printf("ToContinue: %d, %lld, %lld\n", to_continue, heap[0].dist, step2);
#endif
        if (!to_continue) {
            break;
        }
	    if (step == 0) {
			processList(buckets[loIndex][laIndex], lo1e5, la1e5, heap, &len);

	    } else {
			// including both ends
			j = laIndex - step;
			if (j >= 0) {
			    for (i = max(0, loIndex - step); i <= min(loIndex + step, BUCKET_SIZE_LO - 1); i++) {
				    processList(buckets[i][j], lo1e5, la1e5, heap, &len);
			    }	    
			}
			// not including end points
			i = loIndex + step;
			if (i < BUCKET_SIZE_LO) {
			    for (j = max(0, laIndex - step + 1); j <= min(laIndex + step - 1, BUCKET_SIZE_LA - 1); j++) {
				    processList(buckets[i][j], lo1e5, la1e5, heap, &len);
			    }
			}

			// 下面
			j = laIndex + step;
			if (j < BUCKET_SIZE_LA) {
			    for (i = min(BUCKET_SIZE_LO - 1, loIndex + step); i >= max(0, loIndex - step); --i) {
				    processList(buckets[i][j], lo1e5, la1e5, heap, &len);
			    }
			}
			
			// not including end points
			i = loIndex - step;
			if (i >= 0) {
			    for (j = min(BUCKET_SIZE_LA - 1, laIndex + step - 1); j >= min(0, laIndex - step + 1); --j) {
				    processList(buckets[i][j], lo1e5, la1e5, heap, &len);
			    }
			}
	    }

		step++;
	}

	return interpolation(lo, la, heap);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Function to solve the equation a[][] x[] = b[] with the partial-pivoting Gaussian elimination scheme.
 * Copyright (c) Tao Pang 2001.
 */
static void legs (double a[PARAM_NUM][PARAM_NUM], int n, double b1[], double x1[], 
                  double b2[], double x2[], int indx[]) {

    int i,j, idxi, idxj, idxN1;
    double tmp1, tmp2;
    elgs (a, n, indx);

    for(i = 0; i < n-1; ++i) {
        idxi = indx[i];
        for(j = i+1; j < n; ++j) {
            idxj = indx[j];
            b1[idxj] -= a[idxj][i] * b1[idxi];
            b2[idxj] -= a[idxj][i] * b2[idxi];
        }
    }

    idxN1 = indx[n-1];
    x1[n-1] = b1[idxN1] / a[idxN1][n-1];
    x2[n-1] = b2[idxN1] / a[idxN1][n-1];

    for (i = n-2; i>=0; i--) {
        idxi = indx[i];
        
        tmp1 = b1[idxi];
        tmp2 = b2[idxi];

        for (j = i+1; j < n; ++j) {
            tmp1 -= a[idxi][j] * x1[j];
            tmp2 -= a[idxi][j] * x2[j];
        }
        x1[i] = tmp1 / a[idxi][i];
        x2[i] = tmp2 / a[idxi][i];
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Function to perform the partial pivoting Gaussian elimination. a[][] is the original matrix in the input and transformed
 *  matrix plus the pivoting element ratios below the diagonal in the output.  indx[] records the pivoting order.
 *  Copyright (c) Tao Pang 2001.
 */
static void elgs (double a[PARAM_NUM][PARAM_NUM], int n, int indx[]) {
    int i, j, k = 0, itmp;
    double c1, pi, pi1, pj;
    double c[NEIBOR_NUM];
    double tmp;
    int idxi, idxj;

    if (n > PARAM_NUM) {
        printf("The matrix dimension is too large: n = %d, Max: %d \n", n, PARAM_NUM);
        exit(1);
    }

    for (i = 0; i < n; ++i) {
        indx[i] = i;
    }

    // Find the rescaling factors, one from each row
    for (i = 0; i < n; ++i) {
        c1 = 0;
        for (j = 0; j < n; ++j) {
            tmp = fabs(a[i][j]);
            if (tmp > c1) {
                c1 = tmp;
            }
        }
        c[i] = c1;
    }

    // Search the pivoting (largest) element from each column
    for (j = 0; j < n-1; ++j) {
        pi1 = 0;
        for (i = j; i < n; ++i) {
            idxi = indx[i];
            pi = fabs(a[idxi][j])/c[idxi];
            if (pi > pi1) {
                pi1 = pi;
                k = i;
            }
        }

        // Interchange the rows via indx[] to record pivoting order
        itmp = indx[j];
        indx[j] = indx[k];
        indx[k] = itmp;

        idxj = indx[j];
        for (i = j+1; i < n; ++i) {
            idxi = indx[i];
            pj = a[idxi][j]/a[idxj][j];

            // Record pivoting ratios below the diagonal
            a[idxi][j] = pj;

            // Modify other elements accordingly
            for (k = j+1; k < n; ++k) {
                a[idxi][k] = a[idxi][k]-pj * a[idxj][k];
            }
        }
    }
}
static void matrix_mul_aa(double matA[NEIBOR_NUM][PARAM_NUM], 
                          double matAA[PARAM_NUM][PARAM_NUM], int n) {

    int i, j, k;
    double sum;
    // AA = A' * A, n <= NEIBOR_NUM
    for (i = 0; i < PARAM_NUM; i++) {
        for (j = 0; j < PARAM_NUM; j++) {
            sum = 0.0;
            for (k = 0; k < n; k++) {
                sum += matA[k][i] * matA[k][j];
            }
            matAA[i][j] = sum;
        }
    }
}

static void matrix_mul_ab(double matA[NEIBOR_NUM][PARAM_NUM], 
                          const double C[NEIBOR_NUM],  double D[PARAM_NUM], int n) {
    int i, k;
    double sum;
    // D = A' * C
    for (i = 0; i < PARAM_NUM; i++) {
        sum = 0.0;
        for (k = 0; k < n; k++) {
            sum += matA[k][i] * C[k];
        }
        D[i] = sum;
    }
}

static Point interpolation(double lo, double la,
                           Node heap[NEIBOR_NUM]) {

	int i;

	const Point2 *pt2;
	Point result;
	
	int n;

	double matA[NEIBOR_NUM][PARAM_NUM];
	double matAA[PARAM_NUM][PARAM_NUM];
	double los[NEIBOR_NUM], las[NEIBOR_NUM];
	double los1[PARAM_NUM], las1[PARAM_NUM];
	double coeffs1[PARAM_NUM], coeffs2[PARAM_NUM];
	int indx[PARAM_NUM];

	for (i = 0; i < NEIBOR_NUM; i++) {
        pt2 = heap[i].pt;
        //printf("%d %d %d %d\n", pt2->lo0, pt2->la0, pt2->lo1, pt2->la1);
	    matA[i][0] = pt2->lo0;
	    matA[i][1] = pt2->la0;
	    matA[i][2] = 1;

	    los[i] = pt2->lo1;
	    las[i] = pt2->la1;
	}
	n = NEIBOR_NUM;

	matrix_mul_aa(matA, matAA, n);
	matrix_mul_ab(matA, los, los1, n);
	matrix_mul_ab(matA, las, las1, n);
    
    // matAA * coeffs1 = los1
    // matAA * coeffs2 = las1
    legs (matAA, PARAM_NUM, los1, coeffs1, las1, coeffs2, indx);
    result.lo = (coeffs1[0] * lo + coeffs1[1] * la + coeffs1[2]); // * INVERSE_SCALE;
    result.la = (coeffs2[0] * lo + coeffs2[1] * la + coeffs2[2]); // * INVERSE_SCALE;
    
#ifdef PARAM_TEST
    printf("Params: (%9.6f, %9.6f, %15.6f) (%9.6f, %9.6f, %15.6f)\n", coeffs1[0], coeffs1[1], coeffs1[2],
           coeffs2[0], coeffs2[1], coeffs2[2]);
#endif
    
    // printf("After Interpolation: %.5f %.5f\n", result.lo, result.la);
	return result;
}

static int offsets[9] = {0, -900, 900, -1, 1,  -901, -899, 901, 899};
// 将GPS坐标转化成为Sogou坐标
static Point _gpsToSogou(const double lo, const double la) {
	int i;
	Point result;
	
	int n = 9;

	double matA[9][PARAM_NUM];
	double matAA[PARAM_NUM][PARAM_NUM];
	double los[9], las[9];
	double los1[PARAM_NUM], las1[PARAM_NUM];
	double coeffs1[PARAM_NUM], coeffs2[PARAM_NUM];
	int indx[PARAM_NUM];

    // 统一坐标系统
    double lng = (lo - 72.0) * 20;
    double lat = (la - 10.0) * 20;

    int idxLng = (int)lng;
    int idxLat = (int)lat;
    idxLng = min(max(idxLng, 1), 1318);
    idxLat = min(max(idxLat, 1), 898);
    
    int idx = idxLng * 900 + idxLat;

	for (i = 0; i < 9; i++) {
	    int newidx = offsets[i] + idx;

	    int y = newidx / 900; // y为纵向　表示lng
	    int x = newidx % 900; // x为横向, 表示lat
	    PointInt2D* pt = &gps2Sogo[y][x];

	    matA[i][0] = y;
	    matA[i][1] = x;
	    matA[i][2] = 1;

	    los[i] = pt->lo;
	    las[i] = pt->la;
	}


	matrix_mul_aa(matA, matAA, n);
	matrix_mul_ab(matA, los, los1, n);
	matrix_mul_ab(matA, las, las1, n);
    
    // matAA * coeffs1 = los1
    // matAA * coeffs2 = las1
    legs (matAA, PARAM_NUM, los1, coeffs1, las1, coeffs2, indx);
    result.lo = (coeffs1[0] * lng + coeffs1[1] * lat + coeffs1[2]);
    result.la = (coeffs2[0] * lng + coeffs2[1] * lat + coeffs2[2]);
    
#ifdef PARAM_TEST
    printf("Params: (%9.6f, %9.6f, %15.6f) (%9.6f, %9.6f, %15.6f)\n", coeffs1[0], coeffs1[1], coeffs1[2],
           coeffs2[0], coeffs2[1], coeffs2[2]);
#endif
    
    // printf("After Interpolation: %.5f %.5f\n", result.lo, result.la);
	return result;
}

// 将GPS坐标转化成为Google坐标
// GSP->Sogou->Google
// 输入单位为: Degree
static Point _gpsToGoogleMars(const double lo, const double la) {
	Point pt = _gpsToSogou(lo, la);
	pt = _sogouToGoogleMars(pt.lo, pt.la);
	return pt;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ALG_TEST
// Interface to Python
PyObject* geofix_sogouToGoogleMars(PyObject* self, PyObject* args) {
    double lo, la;
    Point pt;

	// 解析参数, 输入为(lo, la)
    if (!PyArg_ParseTuple(args, "dd", &lo, &la)) {
        return NULL;
    }
    
    // printf("input params: %.6f, %.6f\n", lo, la);
	pt = _sogouToGoogleMars(lo, la);
	// printf("input params: %.6f, %.6f\n", pt.lo, pt.la);

    PyObject* tuple = PyTuple_New(2);
    PyTuple_SET_ITEM(tuple, 0, PyFloat_FromDouble(pt.lo * INVERSE_SCALE));
    PyTuple_SET_ITEM(tuple, 1, PyFloat_FromDouble(pt.la * INVERSE_SCALE));
    return tuple;
}

// Interface to Python
PyObject* geofix_gpsToSogou(PyObject* self, PyObject* args) {
    double lo, la;
    Point pt;

	// 解析参数, 输入为(lo, la)
    if (!PyArg_ParseTuple(args, "dd", &lo, &la)) {
        return NULL;
    }

	pt = _gpsToSogou(lo, la);

    PyObject* tuple = PyTuple_New(2);
    PyTuple_SET_ITEM(tuple, 0, PyFloat_FromDouble(pt.lo));
    PyTuple_SET_ITEM(tuple, 1, PyFloat_FromDouble(pt.la));
    return tuple;
}

// Interface to Python
PyObject* geofix_gpsToGoogleMars(PyObject* self, PyObject* args) {
    double lo, la;
    Point pt;

	// 解析参数, 输入为(lo, la)
    if (!PyArg_ParseTuple(args, "dd", &lo, &la)) {
        return NULL;
    }

	pt = _gpsToGoogleMars(lo, la);

    PyObject* tuple = PyTuple_New(2);
    PyTuple_SET_ITEM(tuple, 0, PyFloat_FromDouble(pt.lo * INVERSE_SCALE));
    PyTuple_SET_ITEM(tuple, 1, PyFloat_FromDouble(pt.la * INVERSE_SCALE));
    return tuple;
}


// Interface to Python
PyObject * geofix_sogouLoadData(PyObject *self, PyObject *args) {
    const char *path;
    const char *path2;
    if (!PyArg_ParseTuple(args, "ss", &path, &path2)) {
        return NULL;
    }

	initializeGoogle2Sogou(path);
    initializeGps2Sogou(path2);

    Py_INCREF(Py_None);
    return Py_None;    
}

// Interface to Python
// Our method declararations
static PyMethodDef kGeoFixMethods[] = {
    {"loadAllData", geofix_sogouLoadData, METH_VARARGS,
        "Load Data From Given Path: path1, path2"},
    {"sogouToGoogleMars", geofix_sogouToGoogleMars, METH_VARARGS,
        "sogouToGoogleMars convert from sogou <x, y> to google (lo, la)"},
    {"gpsToSogou", geofix_gpsToSogou, METH_VARARGS,
        "gpsToSogou convert from wgs84(gps) to sogou <x, y>"},
    {"gpsToGoogleMars", geofix_gpsToGoogleMars, METH_VARARGS,
        "gpsToGoogleMars convert from wgs84(gps) <lng, lat> to google (lng, lat)"},        
    {NULL, NULL, 0, NULL}
};

// Interface to Python
// Module initialization
PyMODINIT_FUNC initgeofix_all(void) {
    Py_InitModule("geofix_all", kGeoFixMethods);
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef ALG_TEST 
// Testing Area

void testGps(double lng, double lat, double x, double y) {
	Point pt = _gpsToSogou(lng, lat);
	
	double err1 = pt.lo - x, err2 = pt.la - y;
	double err = sqrt(err1 * err1 + err2 * err2);
#ifndef PARAM_TEST
	printf("| %10.1f | %10.1f | %10.1f | %10.1f | %6d | %6d | %6d | %6c |\n", 
           pt.lo, pt.la, x, y, (int)(err1), (int)(err2 * 10), (int)(err * 10), (err < 2) ? 'Y' : 'N');
#endif
}


void test(double x, double y, double lng, double lat) {
	Point pt = _sogouToGoogleMars(x, y);
	pt.lo = pt.lo / 100000.0;
	pt.la = pt.la / 100000.0;
	
	double err1 = pt.lo - lng, err2 = pt.la - lat;
	double err = sqrt(err1 * err1 + err2 * err2);
#ifndef PARAM_TEST
	printf("| %10.5f | %10.5f | %10.5f | %10.5f | %6d | %6d | %6d | %6c |\n", 
           pt.lo, pt.la, lng, lat, (int)(err1 * 100000), (int)(err2 * 100000), (int)(err * 100000), (err < 0.0001) ? 'Y' : 'N');
#endif
}

int main() {
    initialization("./Sogou2Wgs.txt");
    initializeGps2Sogou("./gps-geomap5.txt");

    printf("-----------------------------------------------------------------------------------\n");
    printf("| %20s | %20s | %15s | %6s | %6s |\n", "Result Geo.    ", "Exp Geo.(lng, lat)  ", "Error(m)   ", "DST(m)", "Pass");
    printf("-----------------------------------------------------------------------------------\n");
    test(12945920,4838036, 116.3, 40);
    test(1.2955019478728717E7, 4854493.340584632,  116.38190182005815,40.113740603544244);
    test(1.2955215777776444E7, 4855249.127650006,  116.3836663456,40.118954345599995);
    
    test(1.2911718286479497E7, 4334918.4958437085, 115.9923643456,36.4342783456);
    test(1.3035199301630903E7, 4295933.870026584,  117.1019833456,36.1509183456);
    test(1.3051735312838482E7, 4194944.964167472,  117.2502993456,35.411171345599996);
    test(1.302790298257175E7,  4019163.729176183,  117.0360653456,34.1069853456);
    test(1.2951339413848642E7, 4845788.30159646,   116.34878,     40.05362); 
    // 清河
    // http://tu.sogou.com:1999/api/smap?lng=116.373011684&lat=40.1188770425&x=12951917.0&y=4843339.0
    test(12951564, 4841884, 116.35080430938172,40.026668567754236);
    test(12951781,4842368, 116.35275641218232,40.03001253753806);
    test(12953170,4803802, 116.36521906042799,39.76315165132224);
    test(12958204,4806723, 116.41048374835975,39.78344148354328);
    
    test(12952357,4842432, 116.35793865632024,40.03046549593787);
    test(12951756,4842408, 116.35253560650165,40.03028827017125);
    
    // 权金城
    test(12951987,4843332, 116.3546098059287,40.03667635288874);
    
    
    printf("-----------------------------------------------------------------------------------\n");
    
    
    testGps(89.100,22.000, 9918674.55299559, 2495374.96035172);
    testGps(110.12,45.34, 1.2258635709840195e7, 5644526.164490255);
    testGps(115.48,45.34, 1.2855314672668554E7, 5644526.164490255);
    // testGps(89.100,22.000, 9918674.55299559, 2495374.96035172);
    
}
#endif
