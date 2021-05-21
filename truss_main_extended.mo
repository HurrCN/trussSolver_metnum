/* 
Numerical Method-01 CCITOnline
Truss Project Solving with Modelica
Developed by Muhammad Hurricane, Faris Ariq Naufal, Laksamana Arria Wibowo

Optimization completed [May 21st 2021]
Problem solved is Moaveni 4e p3/5

A very big thanks to Christopher S. Erwin and Josiah Enrico for helping me
understanding Modelica Language much better than ever!
*/

model TugasProjectTruss_Metnum_Hurricane_3

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
parameter Integer NodesArray[nJoint*2]={ 1, 2,
                                      2, 3,
                                      2, 4,
                                      3, 4};
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
Real u[N];
Real internalForces[N];
//Real normalStress[N];

// =============== MATRICES SETUP ===============
protected
// Properties to find displacement and reaction
Integer boundary[2*size(b,1)]=cat(1,(2*b).-1,2*b);
parameter Integer N=2*nJoint; //define size of global matrices
Real node_i[N]; //
Real node_j[N];
Real G[N,N];
Real g[N,N];
Real G_star[N,N]; //helper, to build global stiffness matrix
Real id[N,N]=identity(N);
Real cos, sin, Length, K[2,2];
// Properties to find internal forces
Real Tc[2,2]; //this is T transpose matrix of {U} = [T]{u}
Real t[N,N];  //this is local global T transpose
Real T[N,N];  //this is global T transpose
Real T_star[N,N];
Real du[N];//Real du;      //this is delta u in order to find the internal forces
Real uix[N];
Real uiy[N];
// float error threshold
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

//Eliminating float error for less than the error value
for i in 1:N loop
  reaction[i]:=if abs(reaction[i])<=error then 0 else reaction[i];
  displacement[i]:=if abs(displacement[i])<=error then 0 else displacement[i];
end for;

// ========================================= INTERNAL FORCES ============================================
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
    //Gain local T transpose matrices
    Tc:=[ cos,   sin;
         -sin,   cos];
    //Transforming the local stiffness to local global T transpose matrices
    t:=zeros(N,N);
    for m in 1:2 loop
      for n in 1:2 loop
        t[2*(Nodes[i,1]-1)+m,2*(Nodes[i,1]-1)+n]:=Tc[m,n];
        t[2*(Nodes[i,2]-1)+m,2*(Nodes[i,2]-1)+n]:=Tc[m,n];
      end for;
    end for;
  
  //accumulate the resulting G to build global T transpose matrix
  T_star:=T+t;
  T:=T_star;
end for;

//Solve the u matrix
u:=T_star*displacement;

for i in 1:nTruss loop
  for j in 1:2 loop
    node_i[j]:=Crd[Nodes[i,1],j];
    node_j[j]:=Crd[Nodes[i,2],j];
  end for;
  Length:=Modelica.Math.Vectors.length(node_j-node_i);
  for k in 1:nJoint loop
    uix[2*k-1]:=u[NodesArray[2*k-1]];
    uiy[2*k]:=u[NodesArray[2*k]];
  du[k]:=(uiy[k]-uix[k]);
  end for;
end for;
internalForces:=(crossArea*young/Length)*du;

for i in 1:N loop
  internalForces[i]:=if abs(internalForces[i])<=error then 0 else internalForces[i];
end for;

end TugasProjectTruss_Metnum_Hurricane_3;
