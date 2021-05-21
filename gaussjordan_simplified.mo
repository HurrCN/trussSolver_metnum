function SubFunction_GaussJordan_Hurricane

//Developed by M. Hurricane May 21st 2021
//This is simplified version code of gauss-jordan calculation

//let's say, we have matrices equation [A]{u} = {b}
input Integer N; //define the order
input Real [N,N] A; //domain matrix
input Real [N] b; //codomain matrix
output Real [N] u; //eigenvector

//define float error
Real error = 10e-10;

algorithm
u:=Modelica.Math.Matrices.solve(A,b); //solving the eigenvector

//Eliminating float error for less than error value
for i in 1:N loop
  if abs(u[i]) <= error then u[i]:= 0;
  end if;
end for;

end SubFunction_GaussJordan_Hurricane;
