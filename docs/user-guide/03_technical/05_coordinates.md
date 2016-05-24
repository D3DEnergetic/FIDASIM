title: Coordinate Systems

#Beam Grid

The coordinates system used depends on the values of the rotation angles (`alpha`, `beta`, `gamma`) and `origin` variables.
If all the rotation angles and `origin` are set to zero then the coordinate system is identical to machine coordinates. 
The angles variables and `origin` define a rotated coordinate system. 
The rotation angles and origin can best be described by an example. 

1. With your right hand point your index finger pointing in the +x direction with your middle finger and thumb pointing in the +y and +z direction respectively.
2. Rotate about your thumb (z-axis) by `alpha` (ccw = +angle, cw = -angle)
3. Rotate about your middle finger (y'-axis) by `beta`
4. Rotate about your index finger (x"-axis) by `gamma`
5. Move your right hand to the `origin`
6. Define `(x|y|z)_(min|max)` by this coordinate system with your index finger being the new +x-axis

