# Computer Graphics

The work in this repository is from a class I took in college. To respect the professor's wishes and to discourage cheating in the professor's future classes, I will not be saying where this class was taken or what the course number was. YOU DO NOT HAVE PERMISSION TO COPY THE CODE IN THIS REPOSITORY.

This class was devoted to learning all about computer graphics. The first half of the class was covered using VTK and the second half focused on OpenGL.

## How to Run the Projects

Clone this repository on your system.

Visit the [VTK website](https://vtk.org/download/) and download the latest version of VTK for your operating system.

Visit the [CMake website](https://cmake.org/download/) and download the latest version of CMake for your operating system.

In each project you will see a src and build folder. The src will contain the source code and a CMakeLists.txt file to help with compiling the project. Depending on the project, the build folder should be empty or contain a vtk data file.

### <u>Mac</u>

Download the latest OS release to ensure OpenGL is up to date.

For VTK projects (projects not labeled OpenGL), modify the CMakeLists.txt file where I have labeled "ADD PATH TO VTK BUILD HERE." You should place the full absolute path to your VTK build in this location. If it's an OpenGL project, ignore this step.

Go to the command line and run _cmake -B/path/to/my/build/folder -S/path/to/my/source/folder_ where the path to the build folder for the project follows the -B and the path to the source folder for the project follows the -S. This should add a Makefile, specific cmake files, and an executable that isn't able to be ran yet to the build folder. It should _not_ modify the src folder at all.

Go into the build folder and run _make_.

Now the executable for the project should be good to run. This can be done by running _./\<project name>.app/Contents/MacOS/\<project name>_ at the command line for VTK projects or _./\<project name>_ for OpenGL projects. Note that most of these projects create images that will appear in the build folder after the executable is ran. The images will have the same name as the project (i.e. the _colors_ project will generate an image called _colors.png_) or if multiple images can be generated they will be named _framexxx.png_ (_xxx_ will be a number).

### <u>Windows</u>

Visit [OpenGL's website](https://www.opengl.org/) and download OpenGL.

Note: I don't currently have a Windows machine so these instructions may be vague/incorrect. Apologies in advance.

For VTK projects (projects not labeled OpenGL), modify the CMakeLists.txt file to accomodate a windows machine. Where I have labeled "ADD PATH TO VTK BUILD HERE" you should add an absolute path to your VTK build folder. For OpenGL projects, the provided CMakeLists.txt file may already work.

Go to the command line and run _cmake -B\path\to\my\build\folder -S\path\to\my\source\folder_ where the path to the build folder for the project follows the -B and the path to the source folder for the project follows the -S. This should add a Makefile, specific cmake files, and an executable that isn't able to be ran yet to the build folder. It should _not_ modify the src folder at all.

Go into the build folder and run _make_.

Now the executable for the project should be good to run. To run the appropriate executable for VTK projects type _.\\\<project name>_ and then tab until you can't anymore (on my Mac this becomes _./\<project name>.app/Contents/MacOS/\<project name>_). For OpenGL projects it will just be _.\\\<project name>_. Note that most of these projects create images that will appear in the build folder after the executable is ran. The images will have the same name as the project (i.e. the _colors_ project will generate an image called _colors.png_) or if multiple images can be generated they will be named _framexxx.png_ (_xxx_ will be a number).
