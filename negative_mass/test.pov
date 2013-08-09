// This is a simple red sphere

// first, the camera position
camera {
  location <2,5,-10>
  look_at <0,0,0>
}

// now, some light
light_source {
  <0,-10,0>
  color rgb <1,1,1>
}

sphere {
  <1,1,1>,5
  pigment {
    marble
    turbulence 1  // full turbulence
    color_map {
      [0.0 color rgb <1,1,0>]
      [0.8 color rgb <1,0,1>]
      [1.0 color rgb <0,1,1>]
    }
  }
}


// the sphere
sphere {
  <0,0,0>, 5
  pigment { color rgb <1,0,0> }
}
