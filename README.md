# Computer-Graphic

I developed a complete 3D software renderer in C++ from scratch without using any GPU APIs. The project focuses on core graphics programming and rendering algorithms.

What I Did:
Wrote a custom rasterizer to draw wireframe and filled triangles using pixel-by-pixel interpolation and depth buffering.

Implemented camera movement and orientation using transformation matrices (glm) for pan, zoom, and orbit controls.

Created an OBJ and MTL file parser to load 3D models and apply material properties.

Designed a lighting system that includes ambient, diffuse, and specular components, calculated per pixel.

Built ray-triangle intersection logic for shadow detection and scene rendering through ray tracing.

Added soft shadows using randomized sampling of an area light source.

Handled input with SDL2 to allow real-time interaction and camera control.

Structured the code to allow switching between multiple rendering modes: wireframe, filled color, shadowed, and ray-traced.

Technologies:
C++, SDL2, GLM, OBJ/MTL parsing, custom rasterization and ray tracing.
