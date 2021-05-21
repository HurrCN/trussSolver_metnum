function SubFunction_GaussJordan_Hurricane

//Developed by M. Hurricane May 21st 2021
//This is simplified version code of gauss-jordan calculation

//let's say, we have matrices equation [A]{u} = {b}
input Real A[8,8]; //domain matrix
input Real b[8]; //codomain matrix
output Real u[8]; //eigenvector

//define float error
protected
Real error = 10e-10;

algorithm
u:=Modelica.Math.Matrices.solve(A,b); //solve to find the eigenvector

//Eliminating float error for less than error value
for i in 1:8 loop
  if abs(u[i]) <= error then u[i]:= 0;
  end if;
end for;

end SubFunction_GaussJordan_Hurricane;
