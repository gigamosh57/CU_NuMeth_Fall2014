################################################
###### Begin Hari Code
################################################

%example1- a script file to solve a system of linear equations
%a script file is a collection of MATLAB commands organized in
%a sequence.
%Note: using ; at the end of a MATLAB command supresses display
%of the result in the command window
%
%This script will not work unless a value for r is provided in the
%command window first
A=[5-r 2+r 1;1 4-r 2;3 2 1-r];
y=[1 2 3]';
%Uses the \ in MATLAB to find inv(A)*y, which is the solution
%to [A].x=y
x=A\y

%example 2- script for solving a quadratic equation
%recall the formula for the roots of a quadratic equation
%we'll evaluate the discriminant first
disc=b^2-4*a*c;
%then the two roots
r1=(-b+sqrt(disc))/2/a;
r2=(-b-sqrt(disc))/2/a;
%then group them in a vector
r=[r1 r2]'


%script example 4.m illustrates the organization of a sequential
%calculation using two functions.
%I want to take a,b and c from elements of the solution
%vector from the function "solution.m", and find the
%roots of the quadratic equation a x^2+b x+c=0
%So, example4.m is a script that manages the above calculation
%by calling 2 functions.  Note that a script maintains complete
%transparency with the command window/global workspace, so
%watch what happens in the workspace window when you run the
%script
%first need to specify a value for the input variable g
g=0.5;
%version 1 - simply call solution and then quad with elements of x
x=solution(g);
r=quad(x(1),x(2),x(3))
%A less efficient version (why do all this when we already have
%the function quad - that helps us see why it is efficient to break
% a problem into steps and use functions for individual steps)
% version 2 - obtain x from solution.m, and r by rewriting whatever
% in quad
    %x=solution(g);
%now define a, b and c from x vector
    %a=x(1);b=x(2);c=x(3);
%solve quadratic equation
    %disc=b^2-4*a*c;
    %r1=(-b+sqrt(disc))/2/a;
    %r2=(-b-sqrt(disc))/2/a;
    %r=[r1 r2]'



% Ex4fun.m
function [r]=ex4fun(g)
%function ex4fun takes a,b and c as elements of the solution
%vector from solution.m and finds the roots of the quadratic
%equation a x^2+b x+c from quad.m.  This is just a function
%version of example4.m
x=solution(g);
a=x(1);b=x(2);c=x(3);
%disc=b^2-4*a*c
%r1=(-b+sqrt(disc))/2/a;
%r2=(-b-sqrt(disc))/2/a;
r=quad(a,b,c);

% Ex4fun2.m
function [r]=ex4fun2(g)
global a c
%function ex4fun2 takes a,b and c as elements of the solution
%vector from the function "solution.m", and finds the
%roots of the quadratic equation a x^2+b x+c using quad2.m
%where a and c are transferred through the global
%we may have thought that a and c must be declared as global
%in the command window first and then specified before this
%function will work.  But as it turns out, because a nd c are defined
%in this function and known to this function, and this function is
%calling quad2, not declaring global in the command window is OK.
%However, note that the global workspace will have no knowledge
%of a, c, or for that matter b and x after this function is called.
x=solution(g);
a=x(1);b=x(2);c=x(3);
%disc=b^2-4*a*c
%r1=(-b+sqrt(disc))/2/a;
%r2=(-b-sqrt(disc))/2/a;
r=quad2(b);

% quad.m
function r=quad(a,b,c)
%same as example2.m, but written as a function
disc=b^2-4*a*c;
%then the two roots
r1=(-b+sqrt(disc))/2/a;
r2=(-b-sqrt(disc))/2/a;
%then group them in a vector
r=[r1 r2]';

% quad2.m
function r=quad2(b)
%modification of quad.m to again illustrate use of global
%no real purpose - just an illustration, showing how to get
%a and c from the command window through global
global a c
disc=b^2-4*a*c;
%then the two roots
r1=(-b+sqrt(disc))/2/a;
r2=(-b-sqrt(disc))/2/a;
%then group them in a vector
r=[r1 r2]';

% solution.m
function [x]=solution(r)
%example1- a function file to solve a system of linear equations,
%same as in example1.m, except that this is a function file rather
%than a script file
%
%  function x=solution(r) will also work, you don't need [...]
%around the output variable from a function if there is only one
%output variable.
%
%the function file must be saved as functionname.m
%To run this function, you need to say something like
%         myx=solution(0.3) in the command window
%where myx is the solution you seek, and 0.3 is your value of r
%
A=[5-r 2+r 1;1 4-r 2;3 2 1-r];
%note that r would have been passed into the function through the list of input
%variables
b=[1 2 3]';
x=A\b;
%see what happens if you use xx=A\b in place of x=A\b
%xx=A\b;
%the output variable identified in the output field of a function
%must be assigned a value in the function
%after this function finishes executing, query the value of b in the
%command window - you will note that the MATLAB workspace or
%command window have no knowledge of variables that are "local"
%to a function

% solution2.m
function [x,detA]=solution2(r)
%example2- a function file to solve a system of linear equations
%and also evaluate the determinant of A.  This function file
%shows how to set up a function file to get more than one variable
%as output.  In this case, you DO NEED a [......] around the list
%of output variables
A=[5-r 2+r 1;1 4-r 2;3 2 1-r];
b=[1 2 3]';
x=A\b;
detA=det(A);
%when you call this function from the command window or a script,
%you will not update both variables unless you use a calling command
%that also has 2 output variables, i.e.
%[myx,mydetA]=solution2(0.5)
%If you use just solution2(0.5) or whatever=solution2(0.5), only the
%first variable that is output from the function will be assigned.
%Look at the workspace window to see what is being updated

% solution3.m
function [x,detA]=solution3(r)
%continuation of previous example to solve a system of linear equations
%but now we want to use b as defined in the command window.  To
%make b from the command window known to the function, declare it
%as a global variable
%global b
A=[5-r 2+r 1;1 4-r 2;3 2 1-r];
%the definition of b has been commented out and is hence missing
%this function cannot access the b entered in the command window
%unless the global declaration is used
%b=[1 2 3]';
x=A\b;
detA=det(A);
%in the command window, you need to first say
%global b 
%this makes b a global variable
%then assign a value for b
%then call the function
%note that if b is not assigned a value that is compatible with what
%the function expects (i.e. a column vector of length 3), the function
%will return with an error
%try b=[1 2 3] and call the function - it will not work
%try b=pi and call the function - it will not work
%try b=[sqrt(2) sqrt(5) sqrt(8)]' - this will work

################################################
###### End Hari Code
################################################

################################################
###### Begin Python Code
################################################
