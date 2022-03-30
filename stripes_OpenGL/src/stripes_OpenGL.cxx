#include <GL/glew.h>    // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "stripes_OpenGL_data.h"

using std::cout;
using std::endl;

#define PHASE3
#define PHASE4
#define PHASE5

unsigned char *GetColorMap(int &textureSize)
{
    unsigned char controlPts[8][3] =
    {
        {  71,  71, 219 },
        {   0,   0,  91 },
        {   0, 255, 255 },
        {   0, 127,   0 },
        { 255, 255,   0 },
        { 255,  96,   0 },
        { 107,   0,   0 },
        { 224,  76,  76 },
    };
    textureSize = 256;
    unsigned char *ptr = new unsigned char[textureSize*3];
    int nControlPts = 8;
    double amountPerPair = ((double)textureSize-1.0)/(nControlPts-1.0);
    for (int i = 0 ; i < textureSize ; i++)
    {
        int lowerControlPt = (int)(i/amountPerPair);
        int upperControlPt = lowerControlPt+1;
        if (upperControlPt >= nControlPts)
            upperControlPt = lowerControlPt; // happens for i == textureSize-1

        double proportion = (i/amountPerPair)-lowerControlPt;
        for (int j = 0 ; j < 3 ; j++)
            ptr[3*i+j] = controlPts[lowerControlPt][j]
                       + proportion*(controlPts[upperControlPt][j]-
                                     controlPts[lowerControlPt][j]);
    }

    return ptr;
}

unsigned char *GetTigerStripes(int &textureSize)
{
    textureSize = 2048;
    unsigned char *ptr = new unsigned char[textureSize];
    int numStripes = 20;
    int valsPerStripe = textureSize/numStripes;
    for (int i = 0 ; i < numStripes ; i++)
    {
        for (int j = 0 ; j < valsPerStripe ; j++)
        {
           int idx = i*valsPerStripe+j;
           if (j < valsPerStripe / 3)
               ptr[idx] = 152;
           else
               ptr[idx] = 255;
        }
    }
    for (int i = numStripes*valsPerStripe ; i < textureSize ; i++)
    {
        ptr[i] = 0;
    }
    return ptr;
}


void _print_shader_info_log(GLuint shader_index) {
  int max_length = 2048;
  int actual_length = 0;
  char shader_log[2048];
  glGetShaderInfoLog(shader_index, max_length, &actual_length, shader_log);
  printf("shader info log for GL index %u:\n%s\n", shader_index, shader_log);
}

GLuint SetupPhase2DataForRendering()
{
  printf("Getting data for Phase 2\n");

  float points[] = {0.5f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    0.0f, 0.5f, 0.0f,
                   -0.5f, 0.0f, 0.0f};

  float colors[] = {1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f,
                    1.0f, 0.0f, 0.0f};

  GLuint indices[] = {0 ,1, 2,
                      1, 2, 3};

  GLuint points_vbo = 0;
  glGenBuffers(1, &points_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(float), points, GL_STATIC_DRAW);

  GLuint colors_vbo = 0;
  glGenBuffers(1, &colors_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(float), colors, GL_STATIC_DRAW);

  GLuint index_vbo;    // Index buffer object
  glGenBuffers( 1, &index_vbo);
  glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_vbo );
  glBufferData( GL_ELEMENT_ARRAY_BUFFER, 12 * sizeof(GLuint), indices, GL_STATIC_DRAW );

  GLuint vao = 0;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_vbo );

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  return vao;
}

const char *phase2VertexShader =
  "#version 400\n"
  "layout (location = 0) in vec3 vertex_position;\n"
  "layout (location = 1) in vec3 vertex_color;\n"
  "out vec3 color;\n"
  "void main() {\n"
  "  color = vertex_color;\n"
  "  gl_Position = vec4(vertex_position, 1.0);\n"
  "}\n";

const char *phase2FragmentShader =
  "#version 400\n"
  "in vec3 color;\n"
  "out vec4 frag_color;\n"
  "void main() {\n"
  "  frag_color = vec4(color, 1.0);\n"
  "}\n";


GLuint SetupPhase345DataForRendering()
{
  printf("Getting data for Phase 3\n");

  // Create the vbo for points
  GLuint p345_points_vbo = 0;
  glGenBuffers(1, &p345_points_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, p345_points_vbo);
  glBufferData(GL_ARRAY_BUFFER, 77535 * sizeof(float), tri_points, GL_STATIC_DRAW);
  
  // Create the vbo for colors
  GLuint p345_colors_vbo = 0;
  glGenBuffers(1, &p345_colors_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, p345_colors_vbo);
  glBufferData(GL_ARRAY_BUFFER, 25845 * sizeof(float), tri_data, GL_STATIC_DRAW);
  
  // Create the vbo for normals
  GLuint p345_normals_vbo = 0;
  glGenBuffers(1, &p345_normals_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, p345_normals_vbo);
  glBufferData(GL_ARRAY_BUFFER, 77535 * sizeof(float), tri_normals, GL_STATIC_DRAW);
  
  /* TEXTURE SECTION */
  
  // Grab info for both textures
  int colorMapSize;
  int tigerStripeSize;
  unsigned char *colorMap = GetColorMap(colorMapSize);
  unsigned char *tigerStripes = GetTigerStripes(tigerStripeSize);
  
  // Create the vbo for textures
  GLuint textures[2];
  // We have 2 textures
  glGenTextures(2, textures);
  
  // Color Map Texture Initialization
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_1D, textures[0]);
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, colorMapSize * sizeof(unsigned char), 0, GL_RGB, GL_UNSIGNED_BYTE, colorMap);
  glGenerateMipmap(GL_TEXTURE_1D);
  
  // Tiger Stripe Texture Initialization
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_1D, textures[1]);
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RED, tigerStripeSize * sizeof(unsigned char), 0, GL_RED, GL_UNSIGNED_BYTE, tigerStripes);
  glGenerateMipmap(GL_TEXTURE_1D);
  
  /* END TEXTURE SECTION*/

  // Create the vbo for indices (use GL_ELEMENT_ARRAY_BUFFER since it's indices)
  GLuint p345_index_vbo;
  glGenBuffers(1, &p345_index_vbo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, p345_index_vbo );
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, 44106 * sizeof(GLuint), tri_indices, GL_STATIC_DRAW );
  
  // Create the vao using the vbos we just created
  GLuint p345_vao = 0;
  glGenVertexArrays(1, &p345_vao);
  glBindVertexArray(p345_vao);
  glBindBuffer(GL_ARRAY_BUFFER, p345_points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, p345_colors_vbo);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, p345_normals_vbo);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, p345_index_vbo );

  
  // Enable the vbos in the vao
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  return p345_vao;
}

const char *phase345VertexShader =
  "#version 400\n"
  "layout (location = 0) in vec3 vertex_position;\n"
  "layout (location = 1) in float vertex_data;\n"
  "layout (location = 2) in vec3 vertex_normal;\n"
  "uniform mat4 MVP;\n"
  "uniform vec3 cameraloc;  // Camera position \n"
  "uniform vec3 lightdir;   // Lighting direction \n"
  "uniform vec4 lightcoeff; // Lighting coeff, Ka, Kd, Ks, alpha\n"
  "out float data;\n"
  "out float shading_amount;\n"
	"out float depth;\n"
  "void main() {\n"
  "  vec4 position = vec4(vertex_position, 1.0);\n"
  "  gl_Position = MVP*position;\n"
  "  data = vertex_data;\n"
	// Calculate depth
	"  depth = gl_Position.z/gl_Position.w;\n"

#ifdef PHASE5
  // Assign shading_amount a value by calculating phong shading
  // camaraloc  : is the location of the camera
  // lightdir   : is the direction of the light
  // lightcoeff : represents a vec4(Ka, Kd, Ks, alpha) from LightingParams of 1F

	// Calculate view direction
	"  vec3 viewDir = normalize(cameraloc - vertex_position);\n"
	// Calculate diffuse value
	"  float LdotN = dot(lightdir, vertex_normal);\n"
	"  float diffuse = max(0, LdotN);\n"
	// Calculate R
	"  vec3 tempN = 2 * LdotN * vertex_normal;\n"
	"  vec3 R = tempN - lightdir;\n"
	// Calculate RdotV and specular value
	"  float RdotV = max(0, dot(R, viewDir));\n"
	"  float specular = pow(RdotV, lightcoeff[3]);\n"
	// Calculate Phong Shading value
	"  shading_amount = lightcoeff[0] + (lightcoeff[1] * diffuse) + (lightcoeff[2] * specular);\n"
#endif

  "}\n";

const char *phase345FragmentShader =
  "#version 400\n"
  "in float data;\n"
  "in float shading_amount;\n"
	"in float depth;\n"
  "out vec4 frag_color;\n"
	// Initialize samplers for textures
	"uniform sampler1D colorMapTexture;\n"
	"uniform sampler1D tigerStripeTexture;\n"
  "void main() {\n"
	// Get the texture coordinates
	"  float texture_data = (data - 1.0)/5.0;\n"
	// Apply tiger stripe texture using depth
	"  vec4 tiger_stripe = texture(tigerStripeTexture, depth);\n"
	// Apply color map texture
	"  frag_color = texture(colorMapTexture, texture_data);\n"
	// Put it all together
	"  frag_color = frag_color * shading_amount * tiger_stripe.x;\n"
  "}\n";

int main() {
  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit()) {
    fprintf(stderr, "ERROR: could not start GLFW3\n");
    return 1;
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow *window = glfwCreateWindow(700, 700, "CIS 441", NULL, NULL);
  if (!window) {
    fprintf(stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    return 1;
  }
  glfwMakeContextCurrent(window);
  // start GLEW extension handler
  glewExperimental = GL_TRUE;
  glewInit();

  // get version info
  const GLubyte *renderer = glGetString(GL_RENDERER); // get renderer string
  const GLubyte *version = glGetString(GL_VERSION);   // version as a string
  printf("Renderer: %s\n", renderer);
  printf("OpenGL version supported %s\n", version);

  // tell GL to only draw onto a pixel if the shape is closer to the viewer
  glEnable(GL_DEPTH_TEST); // enable depth-testing
  glDepthFunc(GL_LESS); // depth-testing interprets a smaller value as "closer"

  GLuint vao = 0;
#ifndef PHASE3
  vao = SetupPhase2DataForRendering();
  const char* vertex_shader = phase2VertexShader;
  const char* fragment_shader = phase2FragmentShader;
#else
  vao = SetupPhase345DataForRendering();
  const char* vertex_shader = phase345VertexShader;
  const char* fragment_shader = phase345FragmentShader;
#endif

  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vs, 1, &vertex_shader, NULL);
  glCompileShader(vs);
  int params = -1;
  glGetShaderiv(vs, GL_COMPILE_STATUS, &params);
  if (GL_TRUE != params) {
    fprintf(stderr, "ERROR: GL shader index %i did not compile\n", vs);
    _print_shader_info_log(vs);
    exit(EXIT_FAILURE);
  }

  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fs, 1, &fragment_shader, NULL);
  glCompileShader(fs);
  glGetShaderiv(fs, GL_COMPILE_STATUS, &params);
  if (GL_TRUE != params) {
    fprintf(stderr, "ERROR: GL shader index %i did not compile\n", fs);
    _print_shader_info_log(fs);
    exit(EXIT_FAILURE);
  }

  GLuint shader_programme = glCreateProgram();
  glAttachShader(shader_programme, fs);
  glAttachShader(shader_programme, vs);
  glLinkProgram(shader_programme);

  glUseProgram(shader_programme);
  
  // Tell OpenGL where our textures live
  GLuint texture1Location = glGetUniformLocation(shader_programme, "colorMapTexture");
  glUniform1i(texture1Location, 0);
  GLuint texture2Location = glGetUniformLocation(shader_programme, "tigerStripeTexture");
  glUniform1i(texture2Location, 1);

#ifdef PHASE3  // Code block for camera transforms
  // Projection matrix : 30° Field of View
  // display size  : 1000x1000
  // display range : 5 unit <-> 200 units
  glm::mat4 Projection = glm::perspective(
      glm::radians(30.0f), (float)1000 / (float)1000,  40.0f, 60.0f);
  glm::vec3 camera(0, 40, 40);
  glm::vec3 origin(0, 0, 0);
  glm::vec3 up(0, 1, 0);
  // Camera matrix
  glm::mat4 View = glm::lookAt(
    camera, // Camera in world space
    origin, // looks at the origin
    up      // and the head is up
  );
  // Model matrix : an identity matrix (model will be at the origin)
  glm::mat4 Model = glm::mat4(1.0f);
  // Our ModelViewProjection : multiplication of our 3 matrices
  glm::mat4 mvp = Projection * View * Model;

  // Get a handle for our "MVP" uniform
  // Only during the initialisation
  GLuint mvploc = glGetUniformLocation(shader_programme, "MVP");
  // Send our transformation to the currently bound shader, in the "MVP" uniform
  // This is done in the main loop since each model will have a different MVP matrix
  // (At least for the M part)
  glUniformMatrix4fv(mvploc, 1, GL_FALSE, &mvp[0][0]);
#endif

#ifdef PHASE5 // Code block for shading parameters
  GLuint camloc = glGetUniformLocation(shader_programme, "cameraloc");
  glUniform3fv(camloc, 1, &camera[0]);
  glm::vec3 lightdir = glm::normalize(camera - origin);   // Direction of light
  GLuint ldirloc = glGetUniformLocation(shader_programme, "lightdir");
  glUniform3fv(ldirloc, 1, &lightdir[0]);
  glm::vec4 lightcoeff(0.3, 0.7, 2.8, 50.5); // Lighting coeff, Ka, Kd, Ks, alpha
  GLuint lcoeloc = glGetUniformLocation(shader_programme, "lightcoeff");
  glUniform4fv(lcoeloc, 1, &lightcoeff[0]);
#endif

  while (!glfwWindowShouldClose(window)) {
    glClearColor(1, 1, 1, 1);
    // wipe the drawing surface clear
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray(vao);
    // Draw triangles

#ifndef PHASE3
    glDrawElements( GL_TRIANGLES, 6, GL_UNSIGNED_INT, NULL );
#else
    // Add correct number of indices
    glDrawElements( GL_TRIANGLES, 44106, GL_UNSIGNED_INT, NULL );
#endif

    // update other events like input handling
    glfwPollEvents();
    // put the stuff we've been drawing onto the display
    glfwSwapBuffers(window);
  }

  // close GL context and any other GLFW resources
  glfwTerminate();
  return 0;
}
