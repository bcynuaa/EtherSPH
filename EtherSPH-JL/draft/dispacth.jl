#=
  @ author: bcynuaa
  @ date: 2023-11-23 22:05:44
  @ description:
 =#

abstract type X end

struct A <: X x end
struct B <: X x end

f(a1::A, a2::A) = a1.x + a2.x;
f(a1::A, b1::B) = a1.x + b1.x + 1;
f(b1::B, a1::A) = b1.x + a1.x + 2;
f(b1::B, b2::B) = b1.x + b2.x + 3;

a1 = A(1);
a2 = A(2);
b1 = B(3);
b2 = B(4);

vec::Vector{X} = [a1, a2, b1, b2];

ma = zeros(4, 4);
for i = 1:4
    for j = 1:4
        ma[i, j] = f(vec[i], vec[j]);
    end
end

@show ma