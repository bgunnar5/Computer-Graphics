#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <cstdlib>

using std::endl;
using std::cerr;
using std::cout;
using std::string;

#include <GL/glew.h>    // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/vec3.hpp>   // glm::vec3
#include <glm/vec4.hpp>   // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>  // glm::translate, glm::rotate, glm::scale

class RenderManager;

void        SetUpDog(int, RenderManager &, int);
void SetUpGround(RenderManager &);
void SetUpSky(RenderManager &);
const char *GetVertexShader();
const char *GetFragmentShader();

// This file is split into four parts:
// - Part 1: code to set up spheres and cylinders
// - Part 2: a "RenderManager" module
// - Part 3: main function
// - Part 4: SetUpDog and the shader programs -- things you modify
//
// It is intended that you will only need to modify code in Part 4.
// That said, you will need functions in Part 2 and should review
// those functions.
// Further, you are encouraged to look through the entire code base.
//


//
//
// PART 1: code to set up spheres and cylinders
//
//

class Triangle
{
  public:
    glm::vec3 v0;
    glm::vec3 v1;
    glm::vec3 v2;
};

std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
{
    std::vector<Triangle> output(4*list.size());
    output.resize(4*list.size());
    for (unsigned int i = 0 ; i < list.size() ; i++)
    {
        Triangle t = list[i];
        glm::vec3 vmid1, vmid2, vmid3;
        vmid1 = (t.v0 + t.v1) / 2.0f;
        vmid2 = (t.v1 + t.v2) / 2.0f;
        vmid3 = (t.v0 + t.v2) / 2.0f;
        output[4*i+0].v0 = t.v0;
        output[4*i+0].v1 = vmid1;
        output[4*i+0].v2 = vmid3;
        output[4*i+1].v0 = t.v1;
        output[4*i+1].v1 = vmid2;
        output[4*i+1].v2 = vmid1;
        output[4*i+2].v0 = t.v2;
        output[4*i+2].v1 = vmid3;
        output[4*i+2].v2 = vmid2;
        output[4*i+3].v0 = vmid1;
        output[4*i+3].v1 = vmid2;
        output[4*i+3].v2 = vmid3;
    }
    return output;
}

void PushVertex(std::vector<float>& coords,
                const glm::vec3& v)
{
  coords.push_back(v.x);
  coords.push_back(v.y);
  coords.push_back(v.z);
}

//
// Sets up a cylinder that is the circle x^2+y^2=1 extruded from
// Z=0 to Z=1.
//
void GetCylinderData(std::vector<float>& coords, std::vector<float>& normals)
{
  int nfacets = 30;
  for (int i = 0 ; i < nfacets ; i++)
  {
    double angle = 3.14159*2.0*i/nfacets;
    double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
    glm::vec3 fnormal(0.0f, 0.0f, 1.0f);
    glm::vec3 bnormal(0.0f, 0.0f, -1.0f);
    glm::vec3 fv0(0.0f, 0.0f, 1.0f);
    glm::vec3 fv1(cos(angle), sin(angle), 1);
    glm::vec3 fv2(cos(nextAngle), sin(nextAngle), 1);
    glm::vec3 bv0(0.0f, 0.0f, 0.0f);
    glm::vec3 bv1(cos(angle), sin(angle), 0);
    glm::vec3 bv2(cos(nextAngle), sin(nextAngle), 0);
    // top and bottom circle vertices
    PushVertex(coords, fv0);
    PushVertex(normals, fnormal);
    PushVertex(coords, fv1);
    PushVertex(normals, fnormal);
    PushVertex(coords, fv2);
    PushVertex(normals, fnormal);
    PushVertex(coords, bv0);
    PushVertex(normals, bnormal);
    PushVertex(coords, bv1);
    PushVertex(normals, bnormal);
    PushVertex(coords, bv2);
    PushVertex(normals, bnormal);
    // curves surface vertices
    glm::vec3 v1normal(cos(angle), sin(angle), 0);
    glm::vec3 v2normal(cos(nextAngle), sin(nextAngle), 0);
    //fv1 fv2 bv1
    PushVertex(coords, fv1);
    PushVertex(normals, v1normal);
    PushVertex(coords, fv2);
    PushVertex(normals, v2normal);
    PushVertex(coords, bv1);
    PushVertex(normals, v1normal);
    //fv2 bv1 bv2
    PushVertex(coords, fv2);
    PushVertex(normals, v2normal);
    PushVertex(coords, bv1);
    PushVertex(normals, v1normal);
    PushVertex(coords, bv2);
    PushVertex(normals, v2normal);
  }
}

//
// Sets up a sphere with equation x^2+y^2+z^2=1
//
void
GetSphereData(std::vector<float>& coords, std::vector<float>& normals)
{
  int recursionLevel = 3;
  std::vector<Triangle> list;
  {
    Triangle t;
    t.v0 = glm::vec3(1.0f,0.0f,0.0f);
    t.v1 = glm::vec3(0.0f,1.0f,0.0f);
    t.v2 = glm::vec3(0.0f,0.0f,1.0f);
    list.push_back(t);
  }
  for (int r = 0 ; r < recursionLevel ; r++)
  {
      list = SplitTriangle(list);
  }

  for (int octant = 0 ; octant < 8 ; octant++)
  {
    glm::mat4 view(1.0f);
    float angle = 90.0f*(octant%4);
    if(angle != 0.0f)
      view = glm::rotate(view, glm::radians(angle), glm::vec3(1, 0, 0));
    if (octant >= 4)
      view = glm::rotate(view, glm::radians(180.0f), glm::vec3(0, 0, 1));
    for(int i = 0; i < list.size(); i++)
    {
      Triangle t = list[i];
      float mag_reci;
      glm::vec3 v0 = view*glm::vec4(t.v0, 1.0f);
      glm::vec3 v1 = view*glm::vec4(t.v1, 1.0f);
      glm::vec3 v2 = view*glm::vec4(t.v2, 1.0f);
      mag_reci = 1.0f / glm::length(v0);
      v0 = glm::vec3(v0.x * mag_reci, v0.y * mag_reci, v0.z * mag_reci);
      mag_reci = 1.0f / glm::length(v1);
      v1 = glm::vec3(v1.x * mag_reci, v1.y * mag_reci, v1.z * mag_reci);
      mag_reci = 1.0f / glm::length(v2);
      v2 = glm::vec3(v2.x * mag_reci, v2.y * mag_reci, v2.z * mag_reci);
      PushVertex(coords, v0);
      PushVertex(coords, v1);
      PushVertex(coords, v2);
      PushVertex(normals, v0);
      PushVertex(normals, v1);
      PushVertex(normals, v2);
    }
  }
}


//
//
// PART 2: RenderManager module
//
//

void _print_shader_info_log(GLuint shader_index) {
  int max_length = 2048;
  int actual_length = 0;
  char shader_log[2048];
  glGetShaderInfoLog(shader_index, max_length, &actual_length, shader_log);
  printf("shader info log for GL index %u:\n%s\n", shader_index, shader_log);
}

class RenderManager
{
  public:
   enum ShapeType
   {
      SPHERE,
      CYLINDER
   };

                 RenderManager();
   void          SetView(glm::vec3 &c, glm::vec3 &, glm::vec3 &);
   void          SetUpGeometry();
   void          SetColor(double r, double g, double b);
   void          Render(ShapeType, glm::mat4 model);
   GLFWwindow   *GetWindow() { return window; };

  private:
   glm::vec3 color;
   GLuint sphereVAO;
   GLuint sphereNumPrimitives;
   GLuint cylinderVAO;
   GLuint cylinderNumPrimitives;
   GLuint mvploc;
   GLuint colorloc;
   GLuint camloc;
   GLuint ldirloc;
   glm::mat4 projection;
   glm::mat4 view;
   GLuint shaderProgram;
   GLFWwindow *window;

   void SetUpWindowAndShaders();
   void MakeModelView(glm::mat4 &);
};

RenderManager::RenderManager()
{
  SetUpWindowAndShaders();
  SetUpGeometry();
  projection = glm::perspective(
        glm::radians(45.0f), (float)1000 / (float)1000,  5.0f, 100.0f);

  // Get a handle for our MVP and color uniforms
  mvploc = glGetUniformLocation(shaderProgram, "MVP");
  colorloc = glGetUniformLocation(shaderProgram, "color");
  camloc = glGetUniformLocation(shaderProgram, "cameraloc");
  ldirloc = glGetUniformLocation(shaderProgram, "lightdir");

  glm::vec4 lightcoeff(0.3, 0.7, 2.8, 50.5); // Lighting coeff, Ka, Kd, Ks, alpha
  GLuint lcoeloc = glGetUniformLocation(shaderProgram, "lightcoeff");
  glUniform4fv(lcoeloc, 1, &lightcoeff[0]);
}

void
RenderManager::SetView(glm::vec3 &camera, glm::vec3 &origin, glm::vec3 &up)
{ 
   glm::mat4 v = glm::lookAt(
                       camera, // Camera in world space
                       origin, // looks at the origin
                       up      // and the head is up
                 );
   view = v; 
   glUniform3fv(camloc, 1, &camera[0]);
   // Direction of light
   glm::vec3 lightdir = glm::normalize(camera - origin);   
   glUniform3fv(ldirloc, 1, &lightdir[0]);
};

void
RenderManager::SetUpWindowAndShaders()
{
  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit()) {
    fprintf(stderr, "ERROR: could not start GLFW3\n");
    exit(EXIT_FAILURE);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  window = glfwCreateWindow(700, 700, "CIS 441", NULL, NULL);
  if (!window) {
    fprintf(stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    exit(EXIT_FAILURE);
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

  const char* vertex_shader = GetVertexShader();
  const char* fragment_shader = GetFragmentShader();

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

  shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, fs);
  glAttachShader(shaderProgram, vs);
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);
}

void RenderManager::SetColor(double r, double g, double b)
{
   color[0] = r;
   color[1] = g;
   color[2] = b;
}

void RenderManager::MakeModelView(glm::mat4 &model)
{
   glm::mat4 modelview = projection * view * model;
   glUniformMatrix4fv(mvploc, 1, GL_FALSE, &modelview[0][0]);
}

void RenderManager::Render(ShapeType st, glm::mat4 model)
{
   int numPrimitives = 0;
   if (st == SPHERE)
   {
      glBindVertexArray(sphereVAO);
      numPrimitives = sphereNumPrimitives;
   }
   else if (st == CYLINDER)
   {
      glBindVertexArray(cylinderVAO);
      numPrimitives = cylinderNumPrimitives;
   }
   MakeModelView(model);
   glUniform3fv(colorloc, 1, &color[0]);
   glDrawElements(GL_TRIANGLES, numPrimitives, GL_UNSIGNED_INT, NULL);
}

void SetUpVBOs(std::vector<float> &coords, std::vector<float> &normals,
               GLuint &points_vbo, GLuint &normals_vbo, GLuint &index_vbo)
{
  int numIndices = coords.size()/3;
  std::vector<GLuint> indices(numIndices);
  for(int i = 0; i < numIndices; i++)
    indices[i] = i;

  points_vbo = 0;
  glGenBuffers(1, &points_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glBufferData(GL_ARRAY_BUFFER, coords.size() * sizeof(float), coords.data(), GL_STATIC_DRAW);

  normals_vbo = 0;
  glGenBuffers(1, &normals_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
  glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);

  index_vbo = 0;    // Index buffer object
  glGenBuffers(1, &index_vbo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);
}

void RenderManager::SetUpGeometry()
{
  std::vector<float> sphereCoords;
  std::vector<float> sphereNormals;
  GetSphereData(sphereCoords, sphereNormals);
  sphereNumPrimitives = sphereCoords.size() / 3;
  GLuint sphere_points_vbo, sphere_normals_vbo, sphere_indices_vbo;
  SetUpVBOs(sphereCoords, sphereNormals, 
            sphere_points_vbo, sphere_normals_vbo, sphere_indices_vbo);

  std::vector<float> cylCoords;
  std::vector<float> cylNormals;
  GetCylinderData(cylCoords, cylNormals);
  cylinderNumPrimitives = cylCoords.size() / 3;
  GLuint cyl_points_vbo, cyl_normals_vbo, cyl_indices_vbo;
  SetUpVBOs(cylCoords, cylNormals, 
            cyl_points_vbo, cyl_normals_vbo, cyl_indices_vbo);

  GLuint vao[2];
  glGenVertexArrays(2, vao);

  glBindVertexArray(vao[SPHERE]);
  sphereVAO = vao[SPHERE];
  glBindBuffer(GL_ARRAY_BUFFER, sphere_points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, sphere_normals_vbo);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_indices_vbo);
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  glBindVertexArray(vao[CYLINDER]);
  cylinderVAO = vao[CYLINDER];
  glBindBuffer(GL_ARRAY_BUFFER, cyl_points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, cyl_normals_vbo);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cyl_indices_vbo);
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
}

//
// PART3: main function
//
void counterDirectionCheck(int min, int max, bool &backwards, const int &counter) {
  // Reverse the counter
  if (counter == max) {
    backwards = true;
  }
  // Keep the counter normal
  else if (counter == min) {
    backwards = false;
  }
}

int main() 
{
  RenderManager rm;
  GLFWwindow *window = rm.GetWindow();

  glm::vec3 origin(0, 0, 0);
  glm::vec3 up(0, 1, 0);
  
  bool backwards = false;
  bool eyeBackwards = false;

  int counter = 135;
  int eyeCounter = -15;
  int cameraPos = 0;
//  int cameraPos = 60;
  while (!glfwWindowShouldClose(window))
  {
    double angle=cameraPos/300.0*2*M_PI;
    cameraPos++;
    
    // Check if we should reverse the counting yet
    counterDirectionCheck(135, 225, backwards, counter);
    
    // Change the way the tail wags
    if (backwards)
      counter--;
    else
    	counter++;
    
    // Check if we should reverse the counting yet
    counterDirectionCheck(15, -15, eyeBackwards, eyeCounter);
    
    // Change the direction of eye movement
    if (eyeBackwards)
      eyeCounter++;
    else
      eyeCounter--;

    glm::vec3 camera(10*sin(angle), 3, 10*cos(angle));
    rm.SetView(camera, origin, up);

    // wipe the drawing surface clear
    glClearColor(0.3, 0.3, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Render all the images
    SetUpDog(counter, rm, eyeCounter);
    SetUpGround(rm);
    SetUpSky(rm);

    // update other events like input handling
    glfwPollEvents();
    // put the stuff we've been drawing onto the display
    glfwSwapBuffers(window);
  }

  // close GL context and any other GLFW resources
  glfwTerminate();
  return 0;
}

glm::mat4 RotateMatrix(float degrees, float x, float y, float z)
{
   glm::mat4 identity(1.0f);
   glm::mat4 rotation = glm::rotate(identity, 
                                    glm::radians(degrees), 
                                    glm::vec3(x, y, z));
   return rotation;
}

glm::mat4 ScaleMatrix(double x, double y, double z)
{
   glm::mat4 identity(1.0f);
   glm::vec3 scale(x, y, z);
   return glm::scale(identity, scale);
}

glm::mat4 TranslateMatrix(double x, double y, double z)
{
   glm::mat4 identity(1.0f);
   glm::vec3 translate(x, y, z);
   return glm::translate(identity, translate);
}

/*
 HEAD SECTION
*/

void SetUpEyeball(glm::mat4 modelSoFar, RenderManager &rm)
{
   // Set up the outer sphere for the eyes
   glm::mat4 scaled10 = ScaleMatrix(0.15, 0.15, 0.15);
   rm.SetColor(1,1,1);
   rm.Render(RenderManager::SPHERE, modelSoFar*scaled10);

   // Set up the inner sphere for the eyes
   glm::mat4 translate = TranslateMatrix(0, 0, 0.95);
   glm::mat4 scaled30 = ScaleMatrix(0.3, 0.3, 0.3);
   rm.SetColor(0,0,0);
   rm.Render(RenderManager::SPHERE, modelSoFar*scaled10*translate*scaled30);
}

void SetUpLip(glm::mat4 modelSoFar, RenderManager &rm) {
  // Set up the scaling and color for each lip
  glm::mat4 scaled15 = ScaleMatrix(0.2, 0.2, 0.2);
  glm::mat4 scaledWide = ScaleMatrix(1.5, 0.4, 1.5);
  rm.SetColor(.8, .78, .62);
  rm.Render(RenderManager::SPHERE, modelSoFar*scaled15*scaledWide);
}

void SetUpTongue(glm::mat4 modelSoFar, RenderManager &rm) {
  // Set up the scaling and color for the tongue
  glm::mat4 scaled15 = ScaleMatrix(.2, .2, .2);
  glm::mat4 scaledWide = ScaleMatrix(1, 0.2, 1.3);
  rm.SetColor(.78, .5, .81);
  rm.Render(RenderManager::SPHERE, modelSoFar*scaled15*scaledWide);
}

void SetUpNose(glm::mat4 modelSoFar, RenderManager &rm) {
  // Set up the scaling and the color for the nose
  glm::mat4 scaled4 = ScaleMatrix(.05, .05, .05);
  rm.SetColor(0, 0, 0);
  rm.Render(RenderManager::SPHERE, modelSoFar*scaled4);
}

void SetUpEar(glm::mat4 modelSoFar, RenderManager &rm) {
  // Set up the scaling and color for each ear
  glm::mat4 scaled80 = ScaleMatrix(.8, .8, .8);
  glm::mat4 scaledWide = ScaleMatrix(.2, .8, .4);
  rm.SetColor(0.403, 0.298, 0.117);
  rm.Render(RenderManager::SPHERE, modelSoFar*scaledWide*scaled80);
}

void SetUpHat(glm::mat4 modelSoFar, RenderManager &rm) {
  // Set up the scaling and color for the bill of the hat
  glm::mat4 scaledBill = ScaleMatrix(.55, .7, .1);
  rm.SetColor(0.709, 0.439, 0.980);
  rm.Render(RenderManager::CYLINDER, modelSoFar*scaledBill);
  
  // Set up the scaling and color for the top of the hat
  glm::mat4 scaledTop = ScaleMatrix(.5, .5, .5);
  glm::mat4 translateTop = TranslateMatrix(0, -.45, -.25);
  rm.SetColor(0.709, 0.439, 0.980);
  rm.Render(RenderManager::SPHERE, modelSoFar*scaledTop*translateTop);
}

void SetUpHead(glm::mat4 modelSoFar, RenderManager &rm, int counter)
{
  // place center of head at X=0, Y=1, Z=.95
  glm::mat4 translate = TranslateMatrix(0, 1, .95);
  
  // Generate the circular skull
  glm::mat4 headTranslate = TranslateMatrix(0, .15, -.45);
  glm::mat4 scaled50 = ScaleMatrix(.5, .5, .5);
  rm.SetColor(.69, .58, .39);
  rm.Render(RenderManager::SPHERE, modelSoFar*translate*headTranslate*scaled50);
  
  // Generate the hat
  glm::mat4 hatTranslate = TranslateMatrix(0, .5, -.45);
  glm::mat4 rotate90 = RotateMatrix(90, 1, 0, 0);
  glm::mat4 tiltHat = RotateMatrix(-30, 1, 0, 0);
  SetUpHat(modelSoFar*translate*hatTranslate*rotate90*tiltHat, rm);
  
  // Generate the left eye
  glm::mat4 leftEyeTranslate = TranslateMatrix(-0.15, 0.25, 0);
  glm::mat4 rotateInFromLeft = RotateMatrix(counter, 1, 1, 1);
  SetUpEyeball(modelSoFar*translate*leftEyeTranslate*rotateInFromLeft, rm);

  // Generate the right eye
  glm::mat4 rightEyeTranslate = TranslateMatrix(0.15, 0.25, 0);
  glm::mat4 rotateInFromRight = RotateMatrix(counter, 1, 1, 1);
  SetUpEyeball(modelSoFar*translate*rightEyeTranslate*rotateInFromRight, rm);

  // Generate the upper lip
  SetUpLip(modelSoFar*translate, rm);
  
  // Generate the bottom lip
  glm::mat4 lowerLipTranslate = TranslateMatrix(0, -.1, 0);
  SetUpLip(modelSoFar*translate*lowerLipTranslate, rm);
  
  // Generate the tongue
  glm::mat4 tongueTranslate = TranslateMatrix(0, -.05, .1);
  SetUpTongue(modelSoFar*translate*tongueTranslate, rm);
  
  // Generate the nose
  glm::mat4 noseTranslate = TranslateMatrix(0, 0.01, .3);
  SetUpNose(modelSoFar*translate*noseTranslate, rm);
  
  // Generate the left ear
  glm::mat4 leftEarTranslate = TranslateMatrix(-.5, -0.2, -0.45);
  SetUpEar(modelSoFar*translate*leftEarTranslate, rm);
  
  // Generate the right ear
  glm::mat4 rightEarTranslate = TranslateMatrix(.5, -0.2, -0.45);
  SetUpEar(modelSoFar*translate*rightEarTranslate, rm);
}

/*
 BODY SECTION
*/

void SetUpLeg(glm::mat4 modelSoFar, RenderManager &rm){
  glm::mat4 scaled50 = ScaleMatrix(.5, .5, .5);
  glm::mat4 scaledWide = ScaleMatrix(.3, .3, 1.4);
  rm.SetColor(.69, .58, .39);
  rm.Render(RenderManager::CYLINDER, modelSoFar*scaled50*scaledWide);
}

void SetUpFoot(glm::mat4 modelSoFar, RenderManager &rm, int white) {
  glm::mat4 scaledFoot = ScaleMatrix(.2, .3, .04);
  if (white)
    rm.SetColor(1, 1, 1);
  else
  	rm.SetColor(0.403, 0.298, 0.117);
  rm.Render(RenderManager::CYLINDER, modelSoFar*scaledFoot);
}

void SetUpNeck(glm::mat4 modelSoFar, RenderManager &rm) {
  glm::mat4 scaled50 = ScaleMatrix(.5, .5, .5);
  glm::mat4 scaledNeck = ScaleMatrix(.3, .3, 1);
  rm.SetColor(.69, .58, .39);
  rm.Render(RenderManager::CYLINDER, modelSoFar*scaled50*scaledNeck);
}

void SetUpTail(glm::mat4 modelSoFar, RenderManager &rm) {
  glm::mat4 scaled20 = ScaleMatrix(.2, .2, .2);
  glm::mat4 scaledTail = ScaleMatrix(.4, .4, 4);
  rm.SetColor(0.403, 0.298, 0.117);
  rm.Render(RenderManager::CYLINDER, modelSoFar*scaled20*scaledTail);
}

void SetUpBody(glm::mat4 modelSoFar, RenderManager &rm, int counter) {
  // place center of body at X=0, Y=0, Z=0
  glm::mat4 translate = TranslateMatrix(0, 0, 0);

  // Generate the body
  glm::mat4 bodyTranslate = TranslateMatrix(0, .15, -.45);
  glm::mat4 scaled80 = ScaleMatrix(.8, .8, .8);
  glm::mat4 scaledWide = ScaleMatrix(.7, .6, 1.5);
  rm.SetColor(.69, .58, .39);
  rm.Render(RenderManager::SPHERE, modelSoFar*translate*bodyTranslate*scaled80*scaledWide);
  
  // Generate the white belly
  glm::mat4 bellyTranslate = TranslateMatrix(0, 0.07, -.45);
  glm::mat4 scaled70 = ScaleMatrix(.7, .7, .7);
  rm.SetColor(1, 1, 1);
  rm.Render(RenderManager::SPHERE, modelSoFar*translate*bellyTranslate*scaled70*scaledWide);
  
  /*
   LEG SECTION
  */
  
  // Generate back left leg
  glm::mat4 backLeftLegTranslate = TranslateMatrix(-.25, .07, -1.2);
  glm::mat4 rotate90 = RotateMatrix(90, 1, 0, 0);
  SetUpLeg(modelSoFar*translate*backLeftLegTranslate*rotate90, rm);
  
  // Generate back right leg
  glm::mat4 backRightLegTranslate = TranslateMatrix(.25, .07, -1.2);
  SetUpLeg(modelSoFar*translate*backRightLegTranslate*rotate90, rm);
  
  // Generate front left leg
  glm::mat4 frontLeftLegTranslate = TranslateMatrix(-.25, .07, .3);
  SetUpLeg(modelSoFar*translate*frontLeftLegTranslate*rotate90, rm);
  
  // Generate front right leg
  glm::mat4 frontRightLegTranslate = TranslateMatrix(.25, .07, .3);
  SetUpLeg(modelSoFar*translate*frontRightLegTranslate*rotate90, rm);
  
  /*
   FOOT SECTION
  */
  
  // Generate back left foot
  glm::mat4 backLeftFootTranslate = TranslateMatrix(-.25, -.6, -1.1);
  SetUpFoot(modelSoFar*translate*backLeftFootTranslate*rotate90, rm, 0);
  
  // Generate back right foot
  glm::mat4 backRightFootTranslate = TranslateMatrix(.25, -.6, -1.1);
  SetUpFoot(modelSoFar*translate*backRightFootTranslate*rotate90, rm, 1);
  
  // Generate back right foot
  glm::mat4 frontLeftFootTranslate = TranslateMatrix(-.25, -.6, .4);
  SetUpFoot(modelSoFar*translate*frontLeftFootTranslate*rotate90, rm, 1);
  
  // Generate back right foot
  glm::mat4 frontRightFootTranslate = TranslateMatrix(.25, -.6, .4);
  SetUpFoot(modelSoFar*translate*frontRightFootTranslate*rotate90, rm, 0);
  
  /*
   TAIL & NECK SECTION
  */
  
  // Generate neck
  glm::mat4 neckTranslate = TranslateMatrix(0, .9, .45);
  glm::mat4 neckRotate = RotateMatrix(35, 1, 0, 0);
  SetUpNeck(modelSoFar*translate*neckTranslate*rotate90*neckRotate, rm);
  
  // Generate tail
  glm::mat4 tailTranslate = TranslateMatrix(0, .3, -1.4);
  glm::mat4 angleTail = RotateMatrix(315, 1, 0, 0);
  glm::mat4 tailRotate = RotateMatrix(counter, 0, 1, 0);
  SetUpTail(modelSoFar*translate*tailTranslate*rotate90*angleTail*tailRotate, rm);
}

/*
 FOOD BOWL SECTION
*/
void SetUpFood(glm::mat4 modelSoFar, RenderManager &rm) {
  glm::mat4 scaledFood = ScaleMatrix(.1, .1, .1);
  float dist;
  glm::mat4 translateFood;
  // The skeleton for the following code can be found at: https://www.geeksforgeeks.org/program-print-circle-pattern/
  // Draw the food in a circular pattern (r represents radius)
  for (int r = 1; r < 4; r++) {
    // Horizontal movement
    for (int i = 0; i <= 2*r; i++) {
      // Vertical movement
      for (int j = 0; j <= 2*r; j++) {
        dist = sqrt((i-r) * (i-r) + (j-r) * (j-r));
        // Dis is within range (radius - 0.5) and (radius + 0.5)
        if (dist > r - 0.5 && dist < r + 0.5) {
          translateFood = TranslateMatrix(i-r, 3-r, j-r);
          rm.SetColor(0.207, 0.172, 0.109);
          rm.Render(RenderManager::SPHERE, modelSoFar*scaledFood*translateFood);
        }
      }
    }
  }
}

void SetUpBowl(glm::mat4 modelSoFar, RenderManager &rm) {
  // place center of bowl at X=0, Y=-.3, Z=1.7
  glm::mat4 translate = TranslateMatrix(0, -.3, 1.7);
  
  glm::mat4 scaled30 = ScaleMatrix(.5, .3, .5);
  glm::mat4 rotate90 = RotateMatrix(90, 1, 0, 0);
  rm.SetColor(1, 1, 1);
  rm.Render(RenderManager::CYLINDER, modelSoFar*translate*scaled30*rotate90);
  
  SetUpFood(modelSoFar*translate, rm);
}

void SetUpGrass(glm::mat4 modelSoFar, RenderManager &rm) {
  // place center of grass at X=0, Y=-.7, Z=0
  glm::mat4 translate = TranslateMatrix(0, -.7, 0);
  
  glm::mat4 scaledGrass = ScaleMatrix(5, .1, 5);
  glm::mat4 rotate90 = RotateMatrix(90, 1, 0, 0);
  rm.SetColor(0, 1, 0);
  rm.Render(RenderManager::CYLINDER, modelSoFar*translate*scaledGrass*rotate90);
}

void SetUpBoneEnds(glm::mat4 modelSoFar, RenderManager &rm, float x, float y, float z) {
  glm::mat4 scaled40 = ScaleMatrix(.13, .13, .13);
  glm::mat4 translate = TranslateMatrix(x, y, z);
  rm.SetColor(0.929, 0.929, 0.929);
  rm.Render(RenderManager::SPHERE, modelSoFar*translate*scaled40);
}

void SetUpBone(glm::mat4 modelSoFar, RenderManager &rm) {
  // place center of bone at X=, Y=, Z=
  glm::mat4 translate = TranslateMatrix(-1.5, -.6, 1);
  
  // Generate the main portion of the bone
  glm::mat4 scaled40 = ScaleMatrix(.1, .1, .4);
  rm.SetColor(0.929, 0.929, 0.929);
  rm.Render(RenderManager::CYLINDER, modelSoFar*translate*scaled40);
  
  // Generate the ends of the bone
  SetUpBoneEnds(modelSoFar*translate, rm, .05, 0, .4);
  SetUpBoneEnds(modelSoFar*translate, rm, -.05, 0, .4);
  SetUpBoneEnds(modelSoFar*translate, rm, .05, 0, -.05);
  SetUpBoneEnds(modelSoFar*translate, rm, -.05, 0, -.05);
}

void SetUpButterflyEnds(glm::mat4 modelSoFar, RenderManager &rm) {
  glm::mat4 scaledEnd = ScaleMatrix(.05, .05, .05);
  rm.SetColor(0, 0, 0);
  rm.Render(RenderManager::SPHERE, modelSoFar*scaledEnd);
}

void SetUpWings(glm::mat4 modelSoFar, RenderManager &rm, string area) {
  glm::mat4 scaledWing;
  if (area == "top")
  	scaledWing = ScaleMatrix(.15, .15, .05);
  else if (area == "bottom")
    scaledWing = ScaleMatrix(.13, .13, .05);
  rm.SetColor(0.984, 0.572, 0.156);
  rm.Render(RenderManager::CYLINDER, modelSoFar*scaledWing);
}

void SetUpButterfly(glm::mat4 modelSoFar, RenderManager &rm) {
  glm::mat4 translate = TranslateMatrix(2, 1, -2);
  
  // Generate the body
  glm::mat4 scaled10 = ScaleMatrix(.05, .05, .4);
  glm::mat4 rotate135 = RotateMatrix(135, 1, 1, 0);
  rm.SetColor(0, 0, 0);
  rm.Render(RenderManager::CYLINDER, modelSoFar*translate*rotate135*scaled10);
  
  // Generate the ends
  SetUpButterflyEnds(modelSoFar*translate*rotate135, rm);
  glm::mat4 bottomTranslate = TranslateMatrix(.2, -.2, -.29);
  SetUpButterflyEnds(modelSoFar*translate*bottomTranslate*rotate135, rm);
  
  // Generate the wings
  glm::mat4 rightTop = TranslateMatrix(.2, 0, 0);
  glm::mat4 rotateNeg150 = RotateMatrix(-150, 1, 0, 0);
  SetUpWings(modelSoFar*translate*rightTop*rotateNeg150, rm, "top");
  glm::mat4 rightBottom = TranslateMatrix(.3, -.1, -.16);
  glm::mat4 rotate70 = RotateMatrix(70, 1, 0, 0);
  SetUpWings(modelSoFar*translate*rightBottom*rotate70, rm, "bottom");
  glm::mat4 leftTop = TranslateMatrix(-0.06, -.03, -.2);
  glm::mat4 rotateNeg125 = RotateMatrix(-125, 1, 0, .5);
  SetUpWings(modelSoFar*translate*leftTop*rotateNeg125, rm, "top");
  glm::mat4 leftBottom = TranslateMatrix(0.08, -.1, -.33);
  glm::mat4 rotate95 = RotateMatrix(95, 1, 0, -.5);
  SetUpWings(modelSoFar*translate*leftBottom*rotate95, rm, "bottom");
}

void
SetUpDog(int counter, RenderManager &rm, int eyeCounter)
{
  glm::mat4 identity(1.0f);

  // Set up all parts of the dog
  SetUpHead(identity, rm, eyeCounter);
  SetUpBody(identity, rm, counter);
}

void SetUpGround(RenderManager &rm) {
  glm::mat4 identity(1.0f);
  
  // Set up the ground and everything on it
  SetUpBowl(identity, rm);
  SetUpGrass(identity, rm);
  SetUpBone(identity, rm);
}

void SetUpSky(RenderManager &rm) {
  glm::mat4 identity(1.0f);
  
  SetUpButterfly(identity, rm);
}
    
const char *GetVertexShader()
{
   static char vertexShader[1024];
   strcpy(vertexShader, 
          "#version 400\n"
          "layout (location = 0) in vec3 vertex_position;\n"
          "layout (location = 1) in vec3 vertex_normal;\n"
          "uniform mat4 MVP;\n"
          "uniform vec3 cameraloc;  // Camera position \n"
          "uniform vec3 lightdir;   // Lighting direction \n"
          "uniform vec4 lightcoeff; // Lighting coeff, Ka, Kd, Ks, alpha\n"
          "out float shading_amount;\n"
          "void main() {\n"
          "  gl_Position = MVP*vec4(vertex_position, 1.0);\n"
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
          "}\n"
         );
   return vertexShader;
}

const char *GetFragmentShader()
{
   static char fragmentShader[1024];
   strcpy(fragmentShader, 
          "#version 400\n"
          "uniform vec3 color;\n"
          "in float shading_amount;\n"
          "out vec4 frag_color;\n"
          "void main() {\n"
          	// Apply shading factor
          	"  frag_color = vec4(color*shading_amount, 1.0);\n"
            "}\n"
         );
   return fragmentShader;
}

