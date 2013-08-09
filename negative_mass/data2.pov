// povray triangle data 
 

#include "colors.inc"
#declare c1 = texture {pigment {color Red}} 
#declare c2 = texture {pigment {color Blue}} 
#declare c3 = texture {pigment {color Green}} 
#declare c4 = texture {pigment {color Yellow}} 
#declare c5 = texture {pigment {color Orange}} 

//light_source {<-4,2,-4> color rgb<0, 1, 1>}
light_source {<0,0,-30> color rgb<1, 1, 1>}
//light_source {<4,2,-4> color White}
//light_source {<1,4,-2> color rgb<1, 0, 1>} 
//light_source {<-4,2,4> color rgb<1, 0, 1>} 

camera { 
 location <4,3,4> 
 look_at <0,2,0>
 rotate <0,110,0>
 } 

sphere { <0,0,0>, 0.1 texture { c2 } }
sphere { <0,0,1>, 0.1 texture { c3 } }
sphere { <0,1,0>, 0.1 texture { c4 } }
sphere { <1,0,0>, 0.1 texture { c5 } }

plane { <0,-10,0>, 0 pigment {checker color White, color Gray}} 

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

rotate <0,0,0> 

translate <0,2,0> 

} 

