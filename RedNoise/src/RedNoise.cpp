#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Colour.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <TextureMap.h>
#include <TexturePoint.h>
#include "glm/vec3.hpp"
#include <ModelTriangle.h>
#include <map>
#include <RayTriangleIntersection.h>
#include <limits>


#define WIDTH 320
#define HEIGHT 240

glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 5.0);
glm::mat3 cameraOrientation = glm::mat3(1.0); // identity matrix //
float focalLength = 2.0;

//This function should return an evenly spaced list of size numberOfValues that contains floating point numbers between from and to.

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    std::vector<float> result;
    float cal = ((to - from) / (numberOfValues - 1));
    float currentValue = from;
    for (int counter = 0; counter < numberOfValues; counter++) {
        result.push_back(currentValue);
        currentValue += cal;
    }
    return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
    std::vector<glm::vec3> result;
    glm::vec3 cal = ((to - from) / float(numberOfValues - 1));
    glm::vec3 currentValue = from;
    for (int counter = 0; counter < numberOfValues; counter++) {
        result.push_back(currentValue);
        currentValue += cal;
    }
    return result;
}


//------------draw Line starts from here -------------//



void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour col, std::vector<std::vector<float>> &depth) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float zDiff = to.depth - from.depth;

    float numberOfValues = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff / numberOfValues;
    float yStepSize = yDiff / numberOfValues;
    float zStepSize = zDiff / numberOfValues;

    for (float i = 0; i <= numberOfValues; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        float z = from.depth + (zStepSize * i);

        uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue);

        if (round(x) >= 0 && round(x) < WIDTH && round(y) >= 0 && round(y) < HEIGHT) {
            if (1 / z < depth[round(y)][round(x)]) {
                window.setPixelColour(round(x), round(y), colour);
                depth[round(y)][round(x)] = 1/z;
            }
        }
    }
}

/*
void drawLineTextureMapping(DrawingWindow &window, CanvasPoint from, CanvasPoint to, TextureMap &textureMap){
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfValues = std::max(abs(xDiff), abs(yDiff));

    float xStepSize = xDiff / numberOfValues;
    float yStepSize = yDiff / numberOfValues;

    float xTexDif = to.texturePoint.x - from.texturePoint.x;
    float yTexDif = to.texturePoint.y - from.texturePoint.y;
    float xTexStepSize = xTexDif / numberOfValues;
    float yTexStepSize = yTexDif / numberOfValues;

    for (float i = 0; i < numberOfValues; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        float xTex = from.texturePoint.x + xTexStepSize * i;
        float yTex = from.texturePoint.y + yTexStepSize * i;

        if (xTex >= 0  && xTex < textureMap.width && yTex >= 0 && yTex < textureMap.height && x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT) {
            uint32_t colour = textureMap.pixels[size_t (xTex) + size_t (yTex) * textureMap.width];
            window.setPixelColour(round(x), round(y), colour);
        }
    }
}
*/

CanvasPoint RandamisedPoint(){ // random
    CanvasPoint canvasPoint = CanvasPoint(rand()%WIDTH, rand()%HEIGHT);
    return canvasPoint;
}
CanvasTriangle RandamisedTriangle() { // random
    CanvasTriangle canvasTriangle(RandamisedPoint(),RandamisedPoint(),RandamisedPoint());
    return canvasTriangle;
}
Colour RandamisedColour() { // random // rand()%256 will give us a random number in the range 0-255
    return Colour(rand() % 255, rand() % 255, rand() % 255);
}



/*
uint32_t cols(Colour colour){
    return  (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);

}

*/


// Drawing the basic line triangle
void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depth) {
    drawLine(window, triangle.v2(), triangle.v1(), colour,depth ), drawLine(window, triangle.v1(), triangle.v0(), colour, depth),drawLine(window, triangle.v0(), triangle.v2(), colour, depth);
}


void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depth){
    if (triangle.v0().y > triangle.v1().y) { std::swap(triangle.vertices[0], triangle.vertices[1]); }
    if (triangle.v0().y > triangle.v2().y) { std::swap(triangle.vertices[0], triangle.vertices[2]); }
    if (triangle.v1().y > triangle.v2().y) { std::swap(triangle.vertices[1], triangle.vertices[2]); }

    CanvasPoint point0 = triangle.vertices[0], point1 = triangle.vertices[1], point2 = triangle.vertices[2];
    float calSlope0 = (point1.x - point0.x) / (point1.y - point0.y);
    float calSlope1 = (point2.x - point1.x) / (point2.y - point1.y);
    float calSlope2 = (point2.x - point0.x) / (point2.y - point0.y);
    float calDepSlope0 = (point1.depth - point0.depth) / (point1.y - point0.y);
    float calDepSlope1 = (point2.depth - point1.depth) / (point2.y - point1.y);
    float calDepSlope2 = (point2.depth - point0.depth) / (point2.y - point0.y);

    for(float y = point0.y; y <= point1.y; y++){
        float xBegin = point0.x + (y - point0.y) * calSlope0;
        float xEnd = point0.x + (y- point0.y) * calSlope2;

        float zBegin = point0.depth + (y - point0.y) * calDepSlope0;
        float zEnd = point0.depth + (y - point0.y) * calDepSlope2;
         xBegin = window.width - xBegin; // for orbit
         xEnd = window.width - xEnd; // for orbit


        drawLine(window, CanvasPoint(xBegin, y , zBegin), CanvasPoint(xEnd, y, zEnd), colour, depth);
    }

    for(float y = point1.y; y <= point2.y; y++) {
        float xBegin = point1.x + (y - point1.y) * calSlope1;
        float xEnd = point0.x + (y - point0.y) * calSlope2;
        float zBegin = point1.depth + (y - point1.y) * calDepSlope1;
        float zEnd = point0.depth + (y - point0.y) * calDepSlope2;
         xBegin = window.width - xBegin; // for orbit
         xEnd = window.width - xEnd;   // for orbit


        drawLine(window, CanvasPoint(xBegin, y, zBegin), CanvasPoint(xEnd, y, zEnd), colour, depth);
    }
}




/*

void textureMapping(DrawingWindow &window, CanvasTriangle triangle, TextureMap &textureMap){
    CanvasPoint firstSideLeft, firstSideRight, secondSideLeft, secondSideRight;
    // sorting the vertices according to y - coordinates------ smallest comes to the bottom and largest comes to the top.

    if (triangle.v0().y > triangle.v1().y) { std::swap(triangle.vertices[0], triangle.vertices[1]); }
    if (triangle.v0().y > triangle.v2().y) { std::swap(triangle.vertices[0], triangle.vertices[2]); }
    if (triangle.v1().y > triangle.v2().y) { std::swap(triangle.vertices[1], triangle.vertices[2]); }

    // loop through all the y - coordinate in this loop
    for (size_t y = triangle.v0().y; y< triangle.v1().y; y++){
        // calculate the interpolation of x coordinate and texture for x, y. So the drawline function can draw the line to fill the triagnle.
        firstSideLeft.x = interpolation(y,triangle.v0().y, triangle.v1().y, triangle.v0().x, triangle.v1().x);
        firstSideLeft.texturePoint.x = interpolation(y,triangle.v0().y, triangle.v1().y, triangle.v0().texturePoint.x, triangle.v1().texturePoint.x);
        firstSideLeft.texturePoint.y = interpolation(y,triangle.v0().y, triangle.v1().y, triangle.v0().texturePoint.y, triangle.v1().texturePoint.y);

        firstSideRight.x = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().x, triangle.v2().x);
        firstSideRight.texturePoint.x = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.x, triangle.v2().texturePoint.x);
        firstSideRight.texturePoint.y = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.y, triangle.v2().texturePoint.y);

        firstSideLeft.y = y; //These are needed. because we need to update y coordinate every iteration.
        firstSideRight.y = y;
        drawLineTextureMapping(window, firstSideRight, firstSideLeft, textureMap);
    }


    for (size_t y = triangle.v1().y; y< triangle.v2().y; y++) {

        secondSideLeft.x = interpolation(y, triangle.v1().y, triangle.v2().y, triangle.v1().x, triangle.v2().x);
        secondSideLeft.texturePoint.x = interpolation(y, triangle.v1().y, triangle.v2().y, triangle.v1().texturePoint.x,triangle.v2().texturePoint.x);
        secondSideLeft.texturePoint.y = interpolation(y, triangle.v1().y, triangle.v2().y, triangle.v1().texturePoint.y,triangle.v2().texturePoint.y);

        secondSideRight.x = interpolation(y, triangle.v0().y, triangle.v2().y, triangle.v0().x, triangle.v2().x);
        secondSideRight.texturePoint.x = interpolation(y, triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.x,triangle.v2().texturePoint.x);
        secondSideRight.texturePoint.y = interpolation(y, triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.y,triangle.v2().texturePoint.y);

        secondSideLeft.y = y;
        secondSideRight.y = y;
        drawLineTextureMapping(window, secondSideRight, secondSideLeft, textureMap);
    }

    drawStrokedTriangle(window, triangle, Colour(255, 255, 255));

}

*/

std::map<std::string, Colour> readMtl(const std::string &filepath) {
    std::ifstream file(filepath);
    std::map<std::string, Colour> mapping;
    std::string lines;
    std::string name;

    while (std::getline(file, lines)) {
        auto line = split(lines, ' ');
        if (line[0] == "newmtl") {name = line[1];}
        else if (line[0] == "Kd") {mapping[name] = Colour(name, stof(line[1]) *255, stof(line[2])*255, stof(line[3])*255);}}
  //  std::cout << "read is done for mtl" << std::endl;
    file.close();
    return mapping;
}



std::vector<ModelTriangle> loadOBJ(const std::string &filepath, const std::map<std::string, Colour> &mapping, float scaling){
    std::vector<glm::vec3> point;
    std::vector<ModelTriangle> modelTriangles;
    std::ifstream file(filepath);
    std::string line;
    Colour col;

    while (std::getline(file, line)){
        std::vector<std::string> element = split(line, ' ');
        if(element[0] == "v"){
            point.push_back(glm::vec3 (std::stof(element[1]) * scaling,std::stof(element[2]) * scaling,std::stof(element[3]) * scaling));
        }else if (element[0] == "f") {
            modelTriangles.push_back(ModelTriangle(point[std::stoi(element[1]) - 1],point[std::stoi(element[2]) - 1],point[std::stoi(element[3]) - 1],col));
        }else if (element[0] == "usemtl"){
            col = mapping.at(element[1]);
        }
    }
   // std::cout << "read is done for obj" << std::endl;
    file.close();
    return modelTriangles;
}



CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float size = 1) {
    glm::vec3  diffPos = cameraOrientation * (vertexPosition - cameraPosition);
    float xDIVz = (diffPos.x / -diffPos.z);
    float yDIVz = (diffPos.y / diffPos.z);
    float u = (focalLength * xDIVz);
    float v = (focalLength * yDIVz);
    CanvasPoint canvaspoint = CanvasPoint(round(u * size + WIDTH/2), round(v * size + HEIGHT/2) );
    canvaspoint.depth = -glm::length(diffPos);
    return canvaspoint;
}

std::vector<std::vector<float>> depthBuffer(int width, int height) {
    std::vector<std::vector<float>> depth;
    std::vector<float> row;
    for(int y=0; y<height; y++) {
        for(int x=0; x<width; x++) {
            row.push_back(INFINITY);
        }
        depth.push_back(row);
    }
    return depth;
}

// Just a wire frame
void drawWireframeView(DrawingWindow &window, std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depth) {
    for(ModelTriangle modelTriangles : modelTriangle) {
        auto v0 = getCanvasIntersectionPoint(modelTriangles.vertices[0], 240);
        auto v1 = getCanvasIntersectionPoint(modelTriangles.vertices[1], 240);
        auto v2 = getCanvasIntersectionPoint(modelTriangles.vertices[2], 240);
        drawStrokedTriangle(window, CanvasTriangle(v0, v1, v2), modelTriangles.colour, depth);
    }
}

// colour with fireframe
void drawWireframeColour(DrawingWindow &window, std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depth) {
    window.clearPixels();
    std::vector<std::vector<float>> dp = depthBuffer(WIDTH, HEIGHT);
    CanvasTriangle canvastriangle;
    for(auto &modelTriangles : modelTriangle) {
        auto v0 = getCanvasIntersectionPoint(modelTriangles.vertices[0], 240);
        auto v1 = getCanvasIntersectionPoint(modelTriangles.vertices[1], 240);
        auto v2 = getCanvasIntersectionPoint(modelTriangles.vertices[2], 240);
        canvastriangle = CanvasTriangle(v0, v1, v2);
        drawFilledTriangle(window, canvastriangle, modelTriangles.colour, dp);
    }
}

void lookAt() {
    glm::vec3 target = glm::vec3(0.0, 0.0, 0.0);
    glm::vec3 camUp = glm::vec3(0.0, 1.0, 0.0);
    glm::vec3 forward = glm::normalize(target - cameraPosition);
    glm::vec3 right = glm::normalize(glm::cross(camUp, forward));
    glm::vec3 up = glm::cross(forward, right);
    cameraOrientation = glm::mat3(right, up, -forward);
}


// add to main two of these for orbit.
//window.clearPixels();
//drawWireframeColour(window, modelTriangles, depth );
//orbit();
void orbit() {
    float orbitSpeed = glm::radians(5.0);
    float cosRotate = cos(orbitSpeed);
    float sinRotate = sin(orbitSpeed);
    glm::mat3 rotationMatrix(
            cosRotate, 0.0, sinRotate,
            0.0,     1.0, 0.0,
            -sinRotate, 0.0, cosRotate
    );
    cameraPosition = rotationMatrix * cameraPosition;
    lookAt();
}

/*

// This is for texture  mapping
void draws(DrawingWindow &window) {
    window.clearPixels();
    TextureMap textureMap("./src/texture.ppm");
    CanvasPoint point0(160, 10);
    CanvasPoint point1(300, 230);
    CanvasPoint point2(10, 150);
    point0.texturePoint = TexturePoint(195, 5);
    point1.texturePoint = TexturePoint(395, 380);
    point2.texturePoint = TexturePoint(65, 330);
        textureMapping(window, CanvasTriangle(point0,point1, point2), textureMap);
}
*/


/* // original draw func
void draw(DrawingWindow &window) {
    window.clearPixels();
    std::vector<float> grayscales = interpolateSingleFloats(255, 0, WIDTH);
    for (size_t y = 0; y < window.height; y++) {
        for (size_t x = 0; x < window.width; x++) {
            //float red = rand() % 256;
            //float green = 0.0;
            //float blue = 0.0;

            float red = grayscales[x];
            float green = grayscales[x];
            float blue = grayscales[x];
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}
*/



void drawRainbow(DrawingWindow &window) {
    window.clearPixels();

    // making the screen to rainbow
    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    std::vector<glm::vec3> beginning = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
    std::vector<glm::vec3> ending = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);

    for (int y = 0; y < HEIGHT; y++) {
        std::vector<glm::vec3> eachrow = interpolateThreeElementValues(beginning[y], ending[y], WIDTH);

        for (int x = 0; x < WIDTH; x++) {
            float red = eachrow[x].r;
            float green = eachrow[x].g;
            float blue = eachrow[x].b;
            uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
            window.setPixelColour(x, y, colour);
        }
    }
}





void handleEvent(SDL_Event event, DrawingWindow &window, std::vector<std::vector<float>> &depth ) {
    float translationSpeed = 0.1;

    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT){
            std::cout << "LEFT" << std::endl;
            cameraPosition.x += translationSpeed;
        }

        else if (event.key.keysym.sym == SDLK_RIGHT){
            std::cout << "RIGHT" << std::endl;
            cameraPosition.x -= translationSpeed;
        }

        else if (event.key.keysym.sym == SDLK_UP){
            std::cout << "UP" << std::endl;
            cameraPosition.y -= translationSpeed;
        }

        else if (event.key.keysym.sym == SDLK_DOWN) {
            std::cout << "DOWN" << std::endl;
            cameraPosition.y += translationSpeed;
        }

        else if (event.key.keysym.sym == SDLK_u) {drawStrokedTriangle(window, RandamisedTriangle(), RandamisedColour(), depth);} // stroked triangle
        else if (event.key.keysym.sym == SDLK_f){ drawFilledTriangle(window, RandamisedTriangle(), RandamisedColour(), depth);} // filled triangle
       // else if (event.key.keysym.sym == SDLK_j){ draws(window);} // texture mapping
        else if(event.key.keysym.sym == SDLK_o){}
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, SDL_FALSE);
    SDL_Event event;

    std::vector<std::vector<float>> depth = depthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> mapping = readMtl("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangles = loadOBJ("./src/cornell-box.obj", mapping, 0.35);

    //drawWireframeView(window, modelTriangles);


     //std::map<std::string, Colour> mtl = readMtl("./src/cornell-box.mtl");
   // std::vector<ModelTriangle> modelTriangle = loadOBJ("./src/cornell-box.obj", 0.35);


    while (true) {




        //drawStrokedTriangle(window, RandamisedTriangle(), RandamisedColour());
        //drawLine(window, CanvasPoint(0, 0), CanvasPoint(WIDTH / 2, HEIGHT / 2), Colour(255, 255, 255));
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window, depth);
        window.clearPixels();
        drawWireframeColour(window, modelTriangles, depth );
        orbit();
        //draw(window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
