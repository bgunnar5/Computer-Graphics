#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;
using std::cout;

double ceil__441(double f)
{
    return ceil(f-0.00001);
}

double floor__441(double f)
{
    return floor(f+0.00001);
}


vtkImageData *NewImage(int height, int width)
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
      double         X[3];
      double         Y[3];
      unsigned char color[3];

  // would some methods for the triangle be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      int width, height;
  // would some methods for accessing and setting pixels be helpful?
};


std::vector<Triangle> GetTriangles(void)
{
   std::vector<Triangle> rv(100);

   unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                  {76,153,0}, {255, 204, 204}, {204, 204, 0}};
   for (int i = 0 ; i < 100 ; i++)
   {
       int idxI = i%10;
       int posI = idxI*100;
       int idxJ = i/10;
       int posJ = idxJ*100;
       int firstPt = (i%3);
       rv[i].X[firstPt] = posI;
       rv[i].Y[firstPt] = posJ;
       rv[i].X[(firstPt+1)%3] = posI;
       rv[i].Y[(firstPt+1)%3] = posJ+99;
       rv[i].X[(firstPt+2)%3] = posI+10*(idxI+1);
       rv[i].Y[(firstPt+2)%3] = posJ+20*(idxJ+1)-50;
       if (i == 5)
           rv[i].Y[firstPt] = -10;
       if (i == 49)
          rv[i].X[(firstPt+2)%3] = posI+20*(idxI+1);
       rv[i].color[0] = colors[i%6][0];
       rv[i].color[1] = colors[i%6][1];
       rv[i].color[2] = colors[i%6][2];
   }

   return rv;
}

void getVertexIndices(Triangle tri, int &right_index, int &bottom_index, int &top_index) {
    /*
    This function finds the index into the X member of the Triangle class for the right vertex,
    bottom left vertex, and the top left vertex.

    ARGS:
        tri - the Triangle object that we're finding indices of the vertices for
        right_index - a passed-by-reference int to store the index of the right vertex
        bottom_index - a passed-by-reference int to store the index of the bottom left vertex
        top_index - a passed-by-reference int to store the index of the top left vertex

    EXAMPLE:
        A triangle with X[3] = {560, 500, 500} and Y[3] = {570, 500, 599} will result in:
        right_index = 0 (corresponding to the point (X[0], Y[0]) or (560, 570) here)
        bottom_index = 1 (corresponding to the point (X[1], Y[1]) or (500, 500) here)
        top_index = 2 (corresponding to the point (X[2], Y[2]) or (500, 599) here)
    */

    // Find the index of the right vertex
   if (tri.X[0] == tri.X[1]) {
       right_index = 2;
   }
   else if (tri.X[1] == tri.X[2]) {
       right_index = 0;
   }
   else if (tri.X[0] == tri.X[2]) {
       right_index = 1;
   }
   
   /*
   Only need to check 2 vertices to find bottom_index and top_index so we figure out which
   2 indices those could be.
   Cases:
        right_index = 0 -> other_index1 = 2, other_index2 = 1
        right_index = 1 -> other_index1 = 0, other_index2 = 2
        right_index = 2 -> other_index1 = 1, other_index2 = 0
   */
   int other_index1 = (right_index+2)%3;
   int other_index2 = (right_index+1)%3;

   if (tri.Y[other_index1] > tri.Y[other_index2]) {
       bottom_index = other_index2;
       top_index = other_index1;
   } else {
       bottom_index = other_index1;
       top_index = other_index2;
   }
}

void initEquationVars(Triangle tri, const int right_index, const int index, double &slope, double &b) {
    /*
    This function initializes the provided slope and b variables so we can use them for future
    equations. We calculate slope using m = rise/run where rise = y2 - y1 and run = x2 - x1. We
    calculate the y-intercept using b = y - mx where we use the right vertex for the x and y
    values.

    ARGS:
        tri - the Triangle object that we're finding a slope and y-intercept for
        right_index - an int representing the index of the right vertex in the triangle's X and 
        Y arrays
        index - an int representing the index of one of the other two vertices besides the
        right vertex (depends on whether we're looking for the top slope or bottom slope)
        slope - a passed-by-reference double to store the value of the slope that we calculate
        b - a passed-by-reference double to store the value of the y-intercept that we calculate

    EXAMPLE:
        Consider a triangle with X[3] = {560, 500, 500} and Y[3] = {570, 500, 599}. The two
        possibilities here are:
        1. Calling initEquationVars(triangle, 0, 1, slope, b) initializes:
            slope = 1.16667
            b = -83.3333
        2. Calling initEquationVars(triangle, 0, 2, slope, b) initializes:
            slope = -0.483333
            b = 840.667
    */

    // Calculate the slope using rise (y2 - y1) over run (x2 - x1)
   double rise = tri.Y[right_index] - tri.Y[index];
   double run = tri.X[right_index] - tri.X[index];
   slope = rise/run;

   // Get the b value using b = y - mx
   b = tri.Y[right_index] - (slope * tri.X[right_index]);
}

int main()
{
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1000*1000;
    for (int i = 0 ; i < npixels*3 ; i++)
        buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;

    // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM

    // Set a variable for the number of triangles we're dealing with
    int num_triangles = 100;

    // Loop through every triangle
    for (int k = 0; k < num_triangles; k++) {
        // Obtain the triangle we're working with
        Triangle tri = triangles[k];

        // Initialize variables to store the index of the 3 vertices
        int right_index;
        int bottom_index;
        int top_index;

        // Assign correct index of right_index, bottom_index, and top_index
        getVertexIndices(tri, right_index, bottom_index, top_index);

        // Initialize variables to be used for the y = mx + b equations
        double top_slope;
        double bottom_slope;
        double top_b;
        double bottom_b;

        // Not dealing with right triangles so we have to split the triangle in half and get
        // two y = mx + b equations (one for top diagonal, another for bottom)
        initEquationVars(tri, right_index, top_index, top_slope, top_b);
        initEquationVars(tri, right_index, bottom_index, bottom_slope, bottom_b);

        // Lambda functions for the two y = mx + b equations
        auto top_y_equation = [top_slope, top_b](double x) {return (top_slope * x) + top_b;};
        auto bottom_y_equation = [bottom_slope, bottom_b](double x) {return (bottom_slope * x) + bottom_b;};

        // Obtain the column min and max values so we know which columns to go over in scanline
        double column_min = ceil__441(tri.X[bottom_index]);
        double column_max = floor__441(tri.X[right_index]);

        // Loop through the correct columns
        for (double i = column_min; i <= column_max; i++) {
            // Check to make sure we're within the scope of the image
            if (i >= screen.width) {
                break;
            }

            // Calculate the bottom and top of the triangle for our specific column
            double bottom_end = ceil__441(bottom_y_equation(i));
            double top_end = floor__441(top_y_equation(i));

            // Loop through the correct rows
            for (double j = bottom_end; j <= top_end; j++) {
                // Check to make sure we're within the scope of the image
                if (j >= screen.height) {
                    break;
                }

                // Calculate the correct index for the pixel
                int index = (int) 3 * (j * screen.width + i);
                if (index < 0) {
                    continue;
                }

                // Color the pixel
                screen.buffer[index + 0] = tri.color[0];
                screen.buffer[index + 1] = tri.color[1];
                screen.buffer[index + 2] = tri.color[2];
            }
        }
    }

    WriteImage(image, "triangles");
}
