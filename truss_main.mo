/* 
Numerical Method-01 CCITOnline
Truss Project with Modelica

Developed by Muhammad Hurricane, Faris Ariq Naufal, Laksamana Arria Wibowo
Problem that is solved is P3/5
*/

model TugasProjectTruss_Metnum_HurricaneFarisLaks

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
                               
//define boundary
parameter Integer a;
equation

end TugasProjectTruss_Metnum_HurricaneFarisLaks;
