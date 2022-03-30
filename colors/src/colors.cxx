#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

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


int main()
{
   std::cerr << "In main!" << endl;

   // Initialize height and width of the image
   int width = 1024;
   int height = 1350;

   // Create the image and get a pointer to the bottom left corner
   vtkImageData *image = NewImage(height, width);
   unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);

   // RGB color value array (index 0 -> color 0; index 1 -> color 127; index 2 -> color 255)
   int colors[] = {0, 127, 255};

   // Initialize the size of each strip and the strip to start at
   int strip_size = 50;
   int strip = 0;

   // Loop through all the rows
   for (int row = 0; row < height; row++) {
      if (row != 0 && row % strip_size == 0){
         strip += 1;
      }

      // Obtain the correct index into the colors array for r, g, and b
      int red_index = strip / 9;
      int green_index = (strip / 3) % 3;
      int blue_index = strip % 3;

      // Set r, g, and b to the correct values based on the index
      int r = colors[red_index];
      int g = colors[green_index];
      int b = colors[blue_index];

      // Loop through all the columns
      for (int col = 0; col < width; col++) {
         // Initialize index variable (3 represents RGB)
         int index = 3 * (row * width + col);

         // Change the color of each pixel appropriately
         buffer[index + 0] = r;
         buffer[index + 1] = g;
         buffer[index + 2] = b;
      }
   }

   WriteImage(image, "colors");
}
