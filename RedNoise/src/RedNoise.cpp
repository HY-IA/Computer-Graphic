#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include "glm/vec3.hpp"

#define WIDTH 320
#define HEIGHT 240

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

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour col) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;

    float numberOfValues = std::max(abs(xDiff), abs(yDiff));
    float xStepSize = xDiff / numberOfValues;
    float yStepSize = yDiff / numberOfValues;

    for (float i = 0; i <= numberOfValues; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        uint32_t colour = (255 << 24) + (int(col.red) << 16) + (int(col.green) << 8) + int(col.blue);
        if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT) {
            window.setPixelColour(round(x), round(y), colour);
        }
    }
}

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



// Drawing the basic line triangle
void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
    drawLine(window, triangle.v2(), triangle.v1(), colour), drawLine(window, triangle.v1(), triangle.v0(), colour),drawLine(window, triangle.v0(), triangle.v2(), colour);
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {


    /*
    CanvasPoint top    = triangle.v0();
    CanvasPoint bottom = triangle.v2();
    CanvasPoint middle = triangle.v1();

    CanvasPoint ver0 = triangle.vertices[0];
    CanvasPoint ver1 = triangle.vertices[1];
    CanvasPoint ver2 = triangle.vertices[2];
*/

    if (triangle.v0().y > triangle.v1().y) { std::swap(triangle.vertices[0], triangle.vertices[1]); }
    if (triangle.v0().y > triangle.v2().y) { std::swap(triangle.vertices[0], triangle.vertices[2]); }
    if (triangle.v1().y > triangle.v2().y) { std::swap(triangle.vertices[1], triangle.vertices[2]); }

    float xPoint0 = triangle.v0().x, yPoint0 = triangle.v0().y;
    float xPoint1 = triangle.v1().x, yPoint1 = triangle.v1().y;
    float xPoint2 = triangle.v2().x, yPoint2 = triangle.v2().y;

    float calSlope0 = (yPoint1 - yPoint0) != 0 ? (xPoint1 - xPoint0) / (yPoint1 - yPoint0) : 0;
    float calSlope1 = (yPoint2 - yPoint1) != 0 ? (xPoint2 - xPoint1) / (yPoint2 - yPoint1) : 0;
    float calSlope2 = (yPoint2 - yPoint0) != 0 ? (xPoint2 - xPoint0) / (yPoint2- yPoint0) : 0;

    for (size_t y = yPoint0; y <= yPoint1; ++y) {
        float xBegin = xPoint0 + (y - yPoint0) * calSlope2;
        float xEnd = xPoint0 + (y - yPoint0) * calSlope0;
        drawLine(window, CanvasPoint(xBegin, y), CanvasPoint(xEnd, y), colour);
    }

    for (size_t y = yPoint1; y <= yPoint2; ++y) {
        float xBegin = xPoint1 + (y - yPoint1) * calSlope1;
        float xEnd = xPoint0 + (y - yPoint0) * calSlope2;
        drawLine(window, CanvasPoint(xBegin, y), CanvasPoint(xEnd, y), colour);
    }

     drawStrokedTriangle(window, triangle, Colour(255, 255, 255));



}





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



void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u) {drawStrokedTriangle(window, RandamisedTriangle(), RandamisedColour());}
        else if (event.key.keysym.sym == SDLK_f){ drawFilledTriangle(window, RandamisedTriangle(), RandamisedColour());}
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, SDL_FALSE);
    SDL_Event event;

    while (true) {
        //drawStrokedTriangle(window, RandamisedTriangle(), RandamisedColour());
        //drawLine(window, CanvasPoint(0, 0), CanvasPoint(WIDTH / 2, HEIGHT / 2), Colour(255, 255, 255));
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window);
        //draw(window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
