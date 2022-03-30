#include <iostream>
#include <string>
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

using std::cerr;
using std::endl;
using std::cout;
using std::min;
using std::max;
using std::string;

double ceil__441(double f)
{
    return ceil(f-0.00001);
}

double floor__441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double X[3];
      double Y[3];
      double Z[3];
      double colors[3][3];

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

std::vector<Triangle> GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("3D_rendering_geometry.vtk");
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
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
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

void sortArrays(double x[], double y[], double z[], double color[][3], int n) {
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
        
        
        // Lerp for the bottom values
        double bottom_y = lerp(rorl.Y[0], rorl.Y[1], bottom_t);
        double bottom_z = lerp(rorl.Z[0], rorl.Z[1], bottom_t);
        double bottom_r = lerp(rorl.colors[0][0], rorl.colors[1][0], bottom_t);
        double bottom_g = lerp(rorl.colors[0][1], rorl.colors[1][1], bottom_t);
        double bottom_b = lerp(rorl.colors[0][2], rorl.colors[1][2], bottom_t);
        
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

            // Check zbuffer for depth and color if appropriate
            if (z > screen.zbuffer[zbuf_index]) {
                screen.buffer[index + 0] = calcColor(r);
                screen.buffer[index + 1] = calcColor(g);
                screen.buffer[index + 2] = calcColor(b);
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

int main()
{
    int width = 1000;
    int height = 1000;
    // Create a 1786 x 1344 image, set the buffer, and fill it with 0s
    vtkImageData *image = NewImage(width, height);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = width*height;
    for (int i = 0 ; i < npixels*3 ; i++)
        buffer[i] = 0;

    double *zbuffer = (double *) malloc(npixels * sizeof(double));
    for (int i = 0; i < npixels; i++)
        zbuffer[i] = -1;

    // Fill triangles vector with appropriate triangles
    std::vector<Triangle> triangles = GetTriangles();

    // Initialize screen appropriately
    Screen screen;
    screen.buffer = buffer;
    screen.width = width;
    screen.height = height;
    screen.zbuffer = zbuffer;


    // Set a variable for the number of triangles we're dealing with
    int num_triangles = triangles.size();

    // Loop through every triangle
    for (int k = 0; k < num_triangles; k++) {
        // Obtain the triangle we're working with
        Triangle tri = triangles[k];

        // Sort the x and y arrays for the triangle in parallel
        sortArrays(tri.X, tri.Y, tri.Z, tri.colors, 3);

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

        // Initialize find variable to keep track of which slope we'll need to calculate in our left/right triangle
        string find;
        // Fill the other two y, z, and rgb values in left/right triangle appropriately
        if (y_value > tri.Y[1]) {
            find = "bottom";
            left.Y[1] = right.Y[1] = tri.Y[1];
            left.Z[1] = right.Z[1] = tri.Z[1];
            left.colors[1][0] = right.colors[1][0] = tri.colors[1][0];
            left.colors[1][1] = right.colors[1][1] = tri.colors[1][1];
            left.colors[1][2] = right.colors[1][2] = tri.colors[1][2];
            left.Y[2] = right.Y[2] = y_value;
            left.Z[2] = right.Z[2] = z;
            left.colors[2][0] = right.colors[2][0] = r;
            left.colors[2][1] = right.colors[2][1] = g;
            left.colors[2][2] = right.colors[2][2] = b_col;
        }
        else {
            find = "top";
            left.Y[1] = right.Y[1] = y_value;
            left.Z[1] = right.Z[1] = z;
            left.colors[1][0] = right.colors[1][0] = r;
            left.colors[1][1] = right.colors[1][1] = g;
            left.colors[1][2] = right.colors[1][2] = b_col;
            left.Y[2] = right.Y[2] = tri.Y[1];
            left.Z[2] = right.Z[2] = tri.Z[1];
            left.colors[2][0] = right.colors[2][0] = tri.colors[1][0];
            left.colors[2][1] = right.colors[2][1] = tri.colors[1][1];
            left.colors[2][2] = right.colors[2][2] = tri.colors[1][2];
        }

        // Rasterize going left triangles
        goingLeft(tri, left, find, slope, b, screen);

        // Rasterize going right triangles
        goingRight(tri, right, find, slope, b, screen);
    }

    // Write to the image
    WriteImage(image, "3D_rendering");
    free(zbuffer);
}
