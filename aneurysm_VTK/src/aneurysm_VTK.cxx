/* Access my movie at: https://www.youtube.com/watch?v=h41GyImcL5I

To create a single image (frame000.png) run the program with #define MOVIE commented out (default).
To create 1000 images run the program with #define MOVIE uncommented. This will create a directory in the current
working directory called imgs and place the 1000 images inside it.
*/
#include <iostream>
#include <string>
#include <cmath>
#include <sys/stat.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#define NORMALS
// #define MOVIE

using std::cerr;
using std::endl;
using std::cout;
using std::min;
using std::max;
using std::string;
using std::to_string;
using std::pow;

/*
IMAGE SECTION
*/

vtkImageData *NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, string filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

/*
CLASS DEFINITION SECTION
*/

class Triangle
{
  public:
      double X[3];
      double Y[3];
      double Z[3];
      double colors[3][3];
      double normals[3][3];
      double shading[3];

  // would some methods for the triangle be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;
      double *zbuffer;
  // would some methods for accessing and setting pixels be helpful?
};

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    Matrix();
    Matrix(double row1[], double row2[], double row3[], double row4[]);
    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

Matrix::Matrix() {};

Matrix::Matrix(double row1[], double row2[], double row3[], double row4[]) {
    // Set values for row 1
    for (int i = 0; i < 4; i++) {
        A[0][i] = row1[i];
    }

    // Set values for row 2
    for (int i = 0; i < 4; i++) {
        A[1][i] = row2[i];
    }

    // Set values for row 3
    for (int i = 0; i < 4; i++) {
        A[2][i] = row3[i];
    }

    // Set values for row 4
    for (int i = 0; i < 4; i++) {
        A[3][i] = row4[i];
    }
}

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "\t(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
} 

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          w[3];
    double          u[3];
    double          v[3];
    double          position[3];
    double          focus[3];
    double          up[3];

    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(int width, int height);
    void            calculateEValues(double cartesian[], double (&arr)[4]);
};

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.8;
         alpha = 50.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;          // The coefficient for ambient lighting
    double Kd;          // The coefficient for diffuse lighting
    double Ks;          // The coefficient for specular lighting
    double alpha;       // The exponent term for specular lighting
};

// LightingParameters lp;

LightingParameters GetLighting(Camera c)
{
    LightingParameters lp;
    lp.lightDir[0] = c.position[0]-c.focus[0];
    lp.lightDir[1] = c.position[1]-c.focus[1];
    lp.lightDir[2] = c.position[2]-c.focus[2];
    double mag = sqrt(lp.lightDir[0]*lp.lightDir[0]
                    + lp.lightDir[1]*lp.lightDir[1]
                    + lp.lightDir[2]*lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}

/*
HELPER FUNCTIONS SECTION
*/

double ceil__441(double f)
{
    return ceil(f-0.00001);
}

double floor__441(double f)
{
    return floor(f+0.00001);
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

std::vector<Triangle> GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("aneurysm_VTK_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

double y_equation(double slope, double x, double b) {
    // y = mx + b
    return (slope * x) + b;
}

void swap(double *xp, double *yp) {
    // Swap two values in an array
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void swapColors(double xp[], double yp[]) {
    // Swap 2D array values
    for (int i = 0; i < 3; i++) {
        double temp = xp[i];
        xp[i] = yp[i];
        yp[i] = temp;
    }
}

void sortArrays(double x[], double y[], double z[], double color[][3], double shading[], int n) {
    // Initialize min index
    int min_idx;

    // Loop through x array until second to last item
    for (int i = 0; i < n-1; i++) {
        // Update min_idx
        min_idx = i;
        // Loop through x array from i+1 to the last item
        for (int j = i+1; j < n; j++) {
            if (x[j] < x[min_idx]) {
                min_idx = j;
            }
        }
        // Sort x and y arrays in parallel
        swap(&x[min_idx], &x[i]);
        swap(&y[min_idx], &y[i]);
        swap(&z[min_idx], &z[i]);
        swap(&shading[min_idx], &shading[i]);
        swapColors(color[min_idx], color[i]);
    }
}

double calcColor(double val) {
    return ceil__441(val*255);
}

double calculateT(double xy, double Axy, double Bxy) {
    return (xy - Axy)/(Bxy - Axy);
}

double lerp(double a, double b, double t) {
    return a + (t * (b - a));
}

void populateOMinFocus(Camera c, double (&arr)[3]) {
    arr[0] = c.position[0] - c.focus[0];
    arr[1] = c.position[1] - c.focus[1];
    arr[2] = c.position[2] - c.focus[2];
}

double normalize(double vector[]) {
    double x_squared = vector[0] * vector[0];
    double y_squared = vector[1] * vector[1];
    double z_squared = vector[2] * vector[2];
    double result = sqrt(x_squared + y_squared + z_squared);
    return result;
}

void crossProduct(double a[], double b[], double (&arr)[3]) {
    // A.y * B.z - A.z * B.y
    arr[0] = (a[1] * b[2]) - (a[2] * b[1]);
    // A.z * B.x - A.x * B.z
    arr[1] = (a[2] * b[0]) - (a[0] * b[2]);
    // A.x * B.y - A.y * B.x
    arr[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

double dotProduct(double a[], double b[]) {
    double x = a[0] * b[0];
    double y = a[1] * b[1];
    double z = a[2] * b[2];
    return x + y + z;
}

void Camera::calculateEValues(double cartesian[], double (&arr)[4]) {
    // Populate e values
    arr[0] = dotProduct(u, cartesian);
    arr[1] = dotProduct(v, cartesian);
    arr[2] = dotProduct(w, cartesian);
}

Matrix Camera::CameraTransform() {
    // Initialize arrays to store the rows
    double row1[4];
    double row2[4];
    double row3[4];
    double row4[4];

    // Initialize cartesian vectors
    double vector1[3] = {1, 0, 0};
    double vector2[3] = {0, 1, 0};
    double vector3[3] = {0, 0, 1};

    // Initialize t vector
    double t[3];
    t[0] = 0 - position[0];
    t[1] = 0 - position[1];
    t[2] = 0 - position[2];

    // Populate the rows
    calculateEValues(vector1, row1);
    calculateEValues(vector2, row2);
    calculateEValues(vector3, row3);
    calculateEValues(t, row4);

    // Fill fourth column
    row1[3] = 0;
    row2[3] = 0;
    row3[3] = 0;
    row4[3] = 1;

    // Create and return the matrix
    Matrix ct(row1, row2, row3, row4);
    return ct;
}

Matrix Camera::ViewTransform() {
    // Initialize arrays to store the rows
    double row1[4];
    double row2[4];
    double row3[4];
    double row4[4];

    // Populate row 1
    row1[0] = 1/tan(angle/2);
    for (int i = 1; i < 4; i++) {
        row1[i] = 0;
    }

    // Populate row 2
    row2[0] = row2[2] = row2[3] = 0;
    row2[1] = 1/tan(angle/2);

    // Populate row 3
    row3[0] = row3[1] = 0;
    row3[2] = (far + near) / (far - near);
    row3[3] = -1;

    // Populate row 4
    row4[0] = row4[1] = row4[3] = 0;
    row4[2] = (2 * far * near) / (far - near);

    // Create and return the matrix
    Matrix vt(row1, row2, row3, row4);
    return vt;
}

Matrix Camera::DeviceTransform(int n, int m) {
    // Initialize arrays to store the rows
    double row1[4];
    double row2[4];
    double row3[4];
    double row4[4];

    // Initialize row 1
    row1[0] = n/2;
    row1[1] = row1[2] = row1[3] = 0;

    // Initialize row 2
    row2[0] = row2[2] = row2[3] = 0;
    row2[1] = m/2;

    // Initialize row 3
    row3[0] = row3[1] = row3[3] = 0;
    row3[2] = 1;

    // Initialize row 4
    row4[0] = n/2;
    row4[1] = m/2;
    row4[2] = 0;
    row4[3] = 1;

    // Create and return the matrix
    Matrix dt(row1, row2, row3, row4);
    return dt;
}

/*
CAMERA FRAME SECTION
*/

void calculateW(Camera &c) {
    // O - focus
    double o_min_focus[3];
    populateOMinFocus(c, o_min_focus);
    // ||O - focus||
    double normalized = normalize(o_min_focus);
    // (O - focus) / ||O - focus||
    c.w[0] = o_min_focus[0] / normalized;
    c.w[1] = o_min_focus[1] / normalized;
    c.w[2] = o_min_focus[2] / normalized;
}

void calculateU(Camera &c) {
    // up x w
    double up_cross_w[3];
    crossProduct(c.up, c.w, up_cross_w);
    // ||up x w||
    double normalized = normalize(up_cross_w);
    // (up x w) / ||up x w||
    c.u[0] = up_cross_w[0] / normalized;
    c.u[1] = up_cross_w[1] / normalized;
    c.u[2] = up_cross_w[2] / normalized;
}

void calculateV(Camera &c) {
    // w x u
    double w_cross_u[3];
    crossProduct(c.w, c.u, c.v);
}

/*
SCANLINE SECTION
*/

void initEquationVars(string type, Triangle tri, double &slope, double &b, double top_y=0, double bottom_y=0) {
    /*
    This function initializes the provided slope and b variables so we can use them for future
    equations. We calculate slope using m = rise/run where rise = y2 - y1 and run = x2 - x1. We
    calculate the y-intercept using b = y - mx where we use the right vertex for the x and y
    values.
    */
    // Initialize variables
    double x;
    double y;
    double rise;
    double run;

    // Assign x and y values appropriately based on which triangle we're dealing with
    if (type == "left") {
        x = tri.X[0];
        y = tri.Y[0];
    }
    else {
        x = tri.X[2];
        y = tri.Y[2];
    }

    // Calculate rise and run based on the type of triangle we're dealing with
    if (type == "full") {
        rise = tri.Y[2] - tri.Y[0];
        run = tri.X[2] - tri.X[0];
    }
    else if (type == "right" && top_y != 0) {
        rise = tri.Y[2] - top_y;
        run = tri.X[2] - tri.X[1];
    }
    else if (type == "right" && top_y == 0) {
        rise = tri.Y[2] - bottom_y;
        run = tri.X[2] - tri.X[1];
    }
    else if (type == "left" && top_y != 0) {
        rise = tri.Y[0] - top_y;
        run = tri.X[0] - tri.X[1];
    }
    else {
        rise = tri.Y[0] - bottom_y;
        run = tri.X[0] - tri.X[1];
    }

    // Calculate the slope using rise (y2 - y1) over run (x2 - x1)
    slope = rise/run;

    // Get the b value using b = y - mx
    b = y - (slope * x);
}

void scanline(Triangle tri, Triangle rorl, Screen screen, double col_min, double col_max, double slope, double b, double other_slope, double other_b, string find) {
    /*
    This function performs the scanline algorithm. Given a column min and max it decides how many columns to 
    loop through. The two slopes and two b values provided allow us to calculate the bottom and top end for each 
    column. The find variable helps designate which y_value calculation we do should be the bottom end and which
    should be the top end.
    */
    // Check to make sure we're within the width of the image
    if (col_max >= screen.width) {
        col_max = screen.width - 1;
    }
    if (col_min < 0) {
        col_min = 0;
    }

    // Loop through the correct columns
    for (double i = col_min; i <= col_max; i++) {
        // Original y corresponds to equation we used in main
        double original_y = y_equation(slope, i, b);
        double other_y = y_equation(other_slope, i, other_b);

        // Initialize bottom_end and top_end variables
        double bottom_end;
        double top_end;

        /* 
        If we needed to calculate the bottom slope, set top_end to be the y_value calculated using the slope and
        b value from main and bottom_end to be the y_value calculated using the slope and b value we got from 
        goingRight/goingLeft.
        */
        if (find == "bottom") {
            bottom_end = ceil__441(other_y);
            top_end = floor__441(original_y);
        }
        /* 
        If we needed to calculate the top slope, set top_end to be the y_value calculated using the slope and
        b value we got from goingRight/goingLeft and bottom_end to be the y_value calculated using the slope 
        and b value from main.
        */
        else if (find == "top") {
            bottom_end = ceil__441(original_y);
            top_end = floor__441(other_y);
        }

        // Check to make sure we're within the height of the image
        if (top_end >= screen.height) {
            top_end = screen.height - 1;
        }
        if (bottom_end < 0) {
            bottom_end = 0;
        }

        // Calculate t values for the top and bottom
        double top_t = calculateT(i, rorl.X[0], rorl.X[2]);
        double bottom_t = calculateT(i, rorl.X[0], rorl.X[1]);

        // Lerp for the top values
        double top_y = lerp(rorl.Y[0], rorl.Y[2], top_t);
        double top_z = lerp(rorl.Z[0], rorl.Z[2], top_t);
        double top_r = lerp(rorl.colors[0][0], rorl.colors[2][0], top_t);
        double top_g = lerp(rorl.colors[0][1], rorl.colors[2][1], top_t);
        double top_b = lerp(rorl.colors[0][2], rorl.colors[2][2], top_t);
        double top_shading = lerp(rorl.shading[0], rorl.shading[2], top_t);
        
        // Lerp for the bottom values
        double bottom_y = lerp(rorl.Y[0], rorl.Y[1], bottom_t);
        double bottom_z = lerp(rorl.Z[0], rorl.Z[1], bottom_t);
        double bottom_r = lerp(rorl.colors[0][0], rorl.colors[1][0], bottom_t);
        double bottom_g = lerp(rorl.colors[0][1], rorl.colors[1][1], bottom_t);
        double bottom_b = lerp(rorl.colors[0][2], rorl.colors[1][2], bottom_t);
        double bottom_shading = lerp(rorl.shading[0], rorl.shading[1], bottom_t);

        // Loop through the correct rows
        for (double j = bottom_end; j <= top_end; j++) {
            // Calculate the correct index for the pixel
            int index = (int) 3 * (j * screen.width + i);
            int zbuf_index = (int) j * screen.width + i;
            if (index < 0) {
                continue;
            }
            // Calculate the t value
            double t = calculateT(j, bottom_y, top_y);

            // Lerp for the rest of the values
            double z = lerp(bottom_z, top_z, t);
            double r = lerp(bottom_r, top_r, t);
            double g = lerp(bottom_g, top_g, t);
            double b = lerp(bottom_b, top_b, t);
            double s = lerp(bottom_shading, top_shading, t);

            // Check zbuffer for depth and color if appropriate
            if (z > screen.zbuffer[zbuf_index]) {
                // Calculate shaded rgb values
                double shaded_r = min((double)1, r * s);
                double shaded_g = min((double)1, g * s);
                double shaded_b = min((double)1, b * s);

                screen.buffer[index + 0] = calcColor(shaded_r);
                screen.buffer[index + 1] = calcColor(shaded_g);
                screen.buffer[index + 2] = calcColor(shaded_b);
                screen.zbuffer[zbuf_index] = z;
            }
        }
    }
}

void goingRight(Triangle tri, Triangle right, string find, double slope, double b, Screen screen) {
    // Determine our column min and max values
    double r_col_min = ceil__441(tri.X[1]);
    double r_col_max = floor__441(tri.X[2]);

    // Check if we should even call scanline
    if (r_col_min > r_col_max) {
        return;
    }

    // Initialize equation variables
    double other_slope;
    double other_b;

    /*
    If we need to find the bottom slope, set the top_y value to 0 and bottom_y value to be the lower y coordinate
    in the middle vertex.
    */
    if (find == "bottom") {
        initEquationVars("right", tri, other_slope, other_b, 0, right.Y[1]);
    }
    // If we need to find the top slope, set the top_y value to be the upper y coordinate in the middle vertex.
    else if (find == "top") {
        initEquationVars("right", tri, other_slope, other_b, right.Y[2]);
    }

    // Call the scanline function with our newfound other_slope and other_b, as well as col min/max
    scanline(tri, right, screen, r_col_min, r_col_max, slope, b, other_slope, other_b, find);
}

void goingLeft(Triangle tri, Triangle left, string find, double slope, double b, Screen screen) {
    // Determine our column min and max values
    double l_col_min = ceil__441(tri.X[0]);
    double l_col_max = floor__441(tri.X[1]);

    // Check if we should even call scanline
    if (l_col_min > l_col_max) {
        return;
    }

    // Initialize equation variables
    double other_slope;
    double other_b;

    /*
    If we need to find the bottom slope, set the top_y value to 0 and bottom_y value to be the lower y coordinate
    in the middle vertex.
    */
    if (find == "bottom") {
        initEquationVars("left", tri, other_slope, other_b, 0, left.Y[1]);
    }
    // If we need to find the top slope, set the top_y value to be the upper y coordinate in the middle vertex.
    else if (find == "top") {
        initEquationVars("left", tri, other_slope, other_b, left.Y[2]);
    }

    // Call the scanline function with our newfound other_slope and other_b, as well as col min/max
    scanline(tri, left, screen, l_col_min, l_col_max, slope, b, other_slope, other_b, find);
}

/*
PHONG SHADING
*/
void calculatePhongShading(LightingParameters &lp, double *viewDir, double normal[3][3], int vertex, Triangle &tri) {
    // LdotN and Diffuse
    double LdotN = dotProduct(lp.lightDir, normal[vertex]);
    double diffuse = max((double)0, LdotN);

    // R and RdotV
    double tempN[3];
    for (int i = 0; i < 3; i++) {
        tempN[i] = 2 * LdotN * normal[vertex][i];
    }
    double R[3];
    R[0] = tempN[0] - lp.lightDir[0];
    R[1] = tempN[1] - lp.lightDir[1];
    R[2] = tempN[2] - lp.lightDir[2];
    double RdotVInitial = dotProduct(R, viewDir);
    double RdotV = max((double)0, RdotVInitial);
    double specular = pow(RdotV, lp.alpha);
    tri.shading[vertex] = lp.Ka + (lp.Kd * diffuse) + (lp.Ks * specular);
}

/*
MAIN SECTION
*/

int main()
{
    // Initialize width and height
    int width = 1000;
    int height = 1000;
    // Fill triangles vector with appropriate triangles
    std::vector<Triangle> triangles = GetTriangles();
    int frames;
    #ifdef MOVIE
    frames = 1000;
    #else
    frames = 1;
    #endif
    int count = 0;
    for (int j = 0; j < frames; j++) {
        // Calculate correct camera frame (0, 250, 500, 750) and create camera
        Camera c = GetCamera(j, 1000);
        LightingParameters lp = GetLighting(c);

        // Make the image
        vtkImageData *image = NewImage(width, height);

        // Create and initialize the buffer (fill with 0s)
        unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
        int npixels = width*height;
        for (int i = 0 ; i < npixels*3 ; i++)
            buffer[i] = 0;

        // Create and initialize zbuffer (fill with -1s)
        double *zbuffer = (double *) malloc(npixels * sizeof(double));
        for (int i = 0; i < npixels; i++)
            zbuffer[i] = -1;

        // Fill triangles vector with appropriate triangles
        // std::vector<Triangle> triangles = GetTriangles();

        // Initialize screen appropriately
        Screen screen;
        screen.buffer = buffer;
        screen.width = width;
        screen.height = height;
        screen.zbuffer = zbuffer;

        // Set up camera frame (stored in Camera object)
        calculateW(c);
        calculateU(c);
        calculateV(c);

        // Create the 3 transform matrices
        Matrix ct = c.CameraTransform();
        Matrix vt = c.ViewTransform();
        Matrix dt = c.DeviceTransform(width, height);

        // Compose the matrices
        Matrix temp = dt.ComposeMatrices(ct, vt);
        Matrix m = temp.ComposeMatrices(temp, dt);

        // Set a variable for the number of triangles we're dealing with
        int num_triangles = triangles.size();

        // // Loop through every triangle
        for (int k = 0; k < num_triangles; k++) {
            // Obtain the triangle we're working with
            // int k = 1;
            Triangle tri = triangles[k];
            
            // Get shading value for each vertex (saved in tri.shading)
            for (int i = 0; i < 3; i++) {
                // View Direction
                double viewDir[3];
                double tempVD[3];
                tempVD[0] = c.position[0] - tri.X[i];  
                tempVD[1] = c.position[1] - tri.Y[i];  
                tempVD[2] = c.position[2] - tri.Z[i];
                double mag = normalize(tempVD);
                if (mag > 0) {
                    viewDir[0] = tempVD[0] / mag;
                    viewDir[1] = tempVD[1] / mag;
                    viewDir[2] = tempVD[2] / mag;
                }
                calculatePhongShading(lp, viewDir, tri.normals, i, tri);
            }


            // Use matrix M to transform each vertex of the triangle
            for (int i = 0; i < 3; i++) {
                double ptIn[4] = {tri.X[i], tri.Y[i], tri.Z[i], 1};
                double *ptOut = new double[4];
                m.TransformPoint(ptIn, ptOut);
                tri.X[i] = ptOut[0]/ptOut[3];
                tri.Y[i] = ptOut[1]/ptOut[3];
                tri.Z[i] = ptOut[2]/ptOut[3];
                delete[] ptOut;
            }

            // Sort the x and y arrays for the triangle in parallel
            sortArrays(tri.X, tri.Y, tri.Z, tri.colors, tri.shading, 3);

            // Initialize the slope and b value between the left and right vertices
            double slope;
            double b;
            string type = "full";
            initEquationVars(type, tri, slope, b);

            // Calculate the point of intersection so we can split the triangle into 2
            double y_value = y_equation(slope, tri.X[1], b);

            // Calculate t value for the intersection point
            double t = calculateT(tri.X[1], tri.X[0], tri.X[2]);

            // Lerp for the z and rgb values at intersection
            double z = lerp(tri.Z[0], tri.Z[2], t);
            double r = lerp(tri.colors[0][0], tri.colors[2][0], t);
            double g = lerp(tri.colors[0][1], tri.colors[2][1], t);
            double b_col = lerp(tri.colors[0][2], tri.colors[2][2], t);
            double s = lerp(tri.shading[0], tri.shading[2], t);

            // Initialize a left and right triangle
            Triangle left;
            Triangle right;

            // Initialize left triangle's x values; left x value, middle, middle
            left.X[0] = tri.X[0];
            left.X[1] = tri.X[1];
            left.X[2] = tri.X[1];
            
            // Initialize one of the left triangle's y, z, and rgb values; other two vertices have to be sorted first since we can't sort by x value
            left.Y[0] = tri.Y[0];
            left.Z[0] = tri.Z[0];
            left.colors[0][0] = tri.colors[0][0];
            left.colors[0][1] = tri.colors[0][1];
            left.colors[0][2] = tri.colors[0][2];
            left.shading[0] = tri.shading[0];

            // Initialize right triangles x values; right x value, middle, middle
            right.X[0] = tri.X[2];
            right.X[1] = tri.X[1];
            right.X[2] = tri.X[1];
        
            // Initialize one of the right triangle's y, z, and rgb values; other two vertices have to be sorted first since we can't sort by x value
            right.Y[0] = tri.Y[2];
            right.Z[0] = tri.Z[2];
            right.colors[0][0] = tri.colors[2][0];
            right.colors[0][1] = tri.colors[2][1];
            right.colors[0][2] = tri.colors[2][2];
            right.shading[0] = tri.shading[2];

            // Initialize find variable to keep track of which slope we'll need to calculate in our left/right triangle
            string find;
            // Fill the other two y, z, and rgb values in left/right triangle appropriately
            if (y_value > tri.Y[1]) {
                find = "bottom";
                left.Y[1] = right.Y[1] = tri.Y[1];
                left.Z[1] = right.Z[1] = tri.Z[1];
                left.shading[1] = right.shading[1] = tri.shading[1];
                left.colors[1][0] = right.colors[1][0] = tri.colors[1][0];
                left.colors[1][1] = right.colors[1][1] = tri.colors[1][1];
                left.colors[1][2] = right.colors[1][2] = tri.colors[1][2];
                left.Y[2] = right.Y[2] = y_value;
                left.Z[2] = right.Z[2] = z;
                left.shading[2] = right.shading[2] = s;
                left.colors[2][0] = right.colors[2][0] = r;
                left.colors[2][1] = right.colors[2][1] = g;
                left.colors[2][2] = right.colors[2][2] = b_col;
            }
            else {
                find = "top";
                left.Y[1] = right.Y[1] = y_value;
                left.Z[1] = right.Z[1] = z;
                left.shading[1] = right.shading[1] = s;
                left.colors[1][0] = right.colors[1][0] = r;
                left.colors[1][1] = right.colors[1][1] = g;
                left.colors[1][2] = right.colors[1][2] = b_col;
                left.Y[2] = right.Y[2] = tri.Y[1];
                left.Z[2] = right.Z[2] = tri.Z[1];
                left.shading[2] = right.shading[2] = tri.shading[1];
                left.colors[2][0] = right.colors[2][0] = tri.colors[1][0];
                left.colors[2][1] = right.colors[2][1] = tri.colors[1][1];
                left.colors[2][2] = right.colors[2][2] = tri.colors[1][2];
            }

            // Rasterize going left triangles
            goingLeft(tri, left, find, slope, b, screen);

            // Rasterize going right triangles
            goingRight(tri, right, find, slope, b, screen);
        }

    // Create name of file
    char str[128];
    #ifdef MOVIE
    if (count == 0) {
        mkdir("imgs", 0777);
    }
    count++;
    sprintf(str, "./imgs/frame%03d", j);
    #else
    sprintf(str, "frame%03d", j);
    #endif

    // Write to the image
    WriteImage(image, str);
    free(zbuffer);
    }
}
