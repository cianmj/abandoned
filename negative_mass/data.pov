// povray triangle data 
 

#include "colors.inc"
#declare c1 = texture {pigment {color Red}} 

light_source {<0,0,0> color White}
light_source {<-10,5,-50> color rgb<1, 1, 1>}
light_source {<0,10,0> color White}
light_source {<-10,-5,50> color rgb<1, 1, 1>} 

camera { 
 location <100,25,-50> 
 look_at <0,0,-50> 
 } 

plane { <0,1,0>, 0 pigment {checker color White, color Gray}} 

mesh { 

triangle{<0,-1,1>,<-1,-1,1>,<-1,0,1> texture { c1 } } 

triangle{<0,-1,1>,<1,-1,1>,<1,0,1> texture { c1 } } 

triangle{<0,-2,0>,<-1,-2,0>,<-1,-1,1> texture { c1 } } 

triangle{<0,-2,0>,<0,-1,1>,<-1,-1,1> texture { c1 } } 

triangle{<0,-2,0>,<1,-2,0>,<1,-1,1> texture { c1 } } 

triangle{<0,-2,0>,<0,-1,1>,<1,-1,1> texture { c1 } } 

triangle{<0,-1,1>,<0,0,-1>,<-1,0,1> texture { c1 } } 

triangle{<0,0,-1>,<0,1,1>,<1,1,1> texture { c1 } } 

triangle{<1,-1,1>,<2,-1,0>,<2,0,0> texture { c1 } } 

triangle{<0,0,-1>,<0,1,1>,<-1,1,1> texture { c1 } } 

triangle{<-1,-1,1>,<-1,0,1>,<-2,0,0> texture { c1 } } 

triangle{<-1,0,1>,<-1,1,1>,<-2,1,0> texture { c1 } } 

triangle{<-1,-2,0>,<-2,-2,0>,<-2,-1,0> texture { c1 } } 

triangle{<0,0,-1>,<-1,0,1>,<-1,1,1> texture { c1 } } 

triangle{<0,1,1>,<0,2,0>,<-1,2,0> texture { c1 } } 

triangle{<0,1,1>,<-1,1,1>,<-1,2,0> texture { c1 } } 

triangle{<-1,-2,0>,<-1,-1,1>,<-2,-1,0> texture { c1 } } 

triangle{<-1,-1,1>,<-2,-1,0>,<-2,0,0> texture { c1 } } 

triangle{<-1,0,1>,<-2,0,0>,<-2,1,0> texture { c1 } } 

triangle{<-1,1,1>,<-1,2,0>,<-2,2,0> texture { c1 } } 

triangle{<-1,1,1>,<-2,1,0>,<-2,2,0> texture { c1 } } 

triangle{<0,-1,1>,<0,0,-1>,<1,0,1> texture { c1 } } 

triangle{<1,-1,1>,<1,0,1>,<2,0,0> texture { c1 } } 

triangle{<1,-2,0>,<2,-2,0>,<2,-1,0> texture { c1 } } 

triangle{<1,0,1>,<1,1,1>,<2,1,0> texture { c1 } } 

triangle{<0,0,-1>,<1,0,1>,<1,1,1> texture { c1 } } 

triangle{<0,1,1>,<0,2,0>,<1,2,0> texture { c1 } } 

triangle{<0,1,1>,<1,1,1>,<1,2,0> texture { c1 } } 

triangle{<1,-2,0>,<1,-1,1>,<2,-1,0> texture { c1 } } 

triangle{<1,0,1>,<2,0,0>,<2,1,0> texture { c1 } } 

triangle{<1,1,1>,<1,2,0>,<2,2,0> texture { c1 } } 

triangle{<1,1,1>,<2,1,0>,<2,2,0> texture { c1 } } 

rotate <-90,0,0>
translate <0,0,0>
} 

