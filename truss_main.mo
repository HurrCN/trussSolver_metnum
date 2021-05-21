/* 
Numerical Method-01 CCITOnline
Truss Project Solving with Modelica

Developed by Muhammad Hurricane, Faris Ariq Naufal, Laksamana Arria Wibowo
Optimization completed [May 21st 2021]
Problem solved is Moaveni 4e p3/5
*/

model TugasProjectTruss_Metnum_Hurricane

// =============== INPUT PROPS & CONDITIONS ===============

//input number of joints given
parameter Integer nJoint=4;
//input number of trusses involved
parameter Integer nTruss=4;
//input cross-sectional area of the trusses, assumed they all are similar
parameter Real crossArea=0.002; //units : m2
//input moduli of elasticity of the trusses
parameter Real young=70e6;      //units : kPa

// =============== DEFINE CONNECTION ===============

//define "Force" as external forces in kilo Newtons
//array format : { U1x, U1y, U2x, U2y, ... , Unx, Uny }  
parameter Real Force[nJoint*2]={ 0, 0,
                                 0, 0,
                                 0, 0,
                                 0, -10};
//define "Nodes" as joint indexing
parameter Integer Nodes[nTruss,2]=[ 1, 2;
                                    2, 3;
                                    2, 4;
                                    3, 4];
//define "Crd" as joint coordinates
parameter Real Crd[nJoint,2]=[ 0,   0;
                               1.5, 0;
                               0,   1;
                               1.5, 1];
                               
// =============== BOUNDARY AND OUTPUTS ===============
  
//define boundary which show joints that will not deflect against the load
parameter Integer b[:]={1,3} ; //means joint 1 and 3 will not deflect

//define output
Real displacement[N];
Real reaction[N];

// =============== MATRICES SETUP ===============
protected
Integer boundary[2*size(b,1)]=cat(1,(2*b).-1,2*b);
parameter Integer N=2*nJoint;
Real node_i[N];
Real node_j[N];
Real G[N,N];
Real g[N,N];
Real G_star[N,N];
Real id[N,N]=identity(N);
Real cos, sin, Length, K[2,2];
Real error=10e-10;

// =============== MAIN ===============
  
algorithm
G:=id;
for i in 1:nTruss loop
  //calling out the nodes' coordinate for each iteration
  for j in 1:2 loop
    node_i[j]:=Crd[Nodes[i,1],j];
    node_j[j]:=Crd[Nodes[i,2],j];
  end for;
  
    //Build up the trigonometric function
    Length:=Modelica.Math.Vectors.length(node_j-node_i); //finding the hypotenuse
    cos:=(node_j[1]-node_i[1])/Length; // cos=adjacent/hypotenuse =sa/mi
    sin:=(node_j[2]-node_i[2])/Length; // sin=opposite/hypotenuse =de/mi
    //Gain local stiffness matrices
    K:=(crossArea*young/Length)*[ cos^2,   cos*sin;
                                  sin*cos,   sin^2];
    //Transforming the local stiffness to local global stiffness matrices
    g:=zeros(N,N);
    for m in 1:2 loop
      for n in 1:2 loop
        g[2*(Nodes[i,1]-1)+m,2*(Nodes[i,1]-1)+n]:=K[m,n];
        g[2*(Nodes[i,2]-1)+m,2*(Nodes[i,2]-1)+n]:=K[m,n];
        g[2*(Nodes[i,2]-1)+m,2*(Nodes[i,1]-1)+n]:=-K[m,n];
        g[2*(Nodes[i,1]-1)+m,2*(Nodes[i,2]-1)+n]:=-K[m,n];
      end for;
    end for;
  
  //accumulate the resulting G to build global stiffness matrix
  G_star:=G+g; 
  G:=G_star;
end for;

//Implementing boundary
for i in boundary loop
  for j in 1:N loop
    G[i,j]:=id[i,j];
  end for;
end for;

//Solving the result based on the formulas given in Moaveni CH 3
displacement:=SubFunction_GaussJordan_Hurricane(G,Force); //taking gaussjordan function to solve the matrices
reaction:=(G_star*displacement)-Force; // R = [K]{U} - {F}

//Eliminating float error for less then the error value
for i in 1:N loop
  reaction[i]:=if abs(reaction[i])<=error then 0 else reaction[i];
  displacement[i]:=if abs(displacement[i])<=error then 0 else displacement[i];
end for;


end TugasProjectTruss_Metnum_Hurricane;
