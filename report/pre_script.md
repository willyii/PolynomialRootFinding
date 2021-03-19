## Page 1 

Hello professor, thank you for participate my defense committee. My name is
Xinlong Yi, My advisor is professor Schroeder.  I'm here for present my project.
The topic of this project is Real-root isolation of polynomials. It's a very
common problem in mathmatics. 

## Page 2

Here is outline today, I will introduce the overview of this project. Then I
will introduce you the theories and algorithms that applied to this project.
Thirdly, how we implemented this project and what we tried to optimize it will
be introduced. Then results and analysis will we introduced. And last, I will
talk about conclusions.

## Page 3

The first is overview, the main purpose of this project is to implement
real-root isolation program based on Budan's theorem and continued fraction.
Both methods can work only on the square-free polynomials. A polynomial is
square-free means it has no repeated root. Therefore, we need use square-free
decomposition before use these two methods.

## Page 4

Then let's go to the theory of this project. This part should start from the
square-free decomposition. We use Yun's algorithm to perform square-free
decomposition. Assume a polynomial is the multiplication of several square-free
polynomials. It forms like P equals to P1 time P2's square times P3's cubic up
to Pk's power of k. Yun's algorithm can return the Pi and its corresponding
repeat times. 

It's noticeable that, this process is based on the success of computation of
Greatest Common Divisor, which we called GCD. This process is very simllar to
the computation of GCD of two integers. The output of previous step is the input
of next step until the result becomes zero. However, since some floating point
computation is inexact in computer, this process is very easy to get error. We
apply interval arithmetic to handle these errors. I will talk about this in
detail later.

## Page 5

After we get square-free polynomials, we could apply Budan's theorem and
continued fraction. Before I begin to introduce these two methods, I'd like to
introduce the Descartes' rule of signs. I talk about this first since both
methods are based on this rule. 

Let's first look at what is sign variation. It defines like given a sequence of
numbers. If a non-zero number have different sign with next non-zero number, we
say this is one sign change or one sign variation. For example, we have 1, 0,
-1000. There is only one sign change of it, since 1 and -1000 have different
sign and zero should be ignored in this sequence.

Then based on the sign variation, Descartes's rule of signs states that given a
polynomial P, the number of sign variation of coefficient sequence minus the
number of positive roots is a nonnegative even integer. Suppose we have a
polynomial whose sign variation of coefficient is 3. Then the possible number of
positive roots are 3 and 1. And if we want negative roots, we change change the
x in P to -x.

## Page 6

Then, based on Descartes' rule of signs, Budan proposed his own statement. Which
can help us bounding the number of real roots in an interval. It states that,
the sign variation of P(x+l) - sign variation of P(x+r) - the number of roots in half open
interval  l to r is a nonnegative even integer. 

Let's look at this example. This polynomial has three roots, 1, 3, 5. The sign
variation of P(x+1) is 2 and sign variation of P(x+4) is 1. Therefore the
possible number of roots in interval 1 to 4 is 1. Actually, it refers to 3. And
since 1 on open end of this interval, it doesn't count.

## Page 7

Then we combine Budan's theorem with bisection method. At first, we get the
upper bound of roots and use -upperbound as lower bound of it. Then we cut the
interval half by half. 

There are two conditions that can terminate this process. The first one is there is
only one real root in searching interval, we isolate it success. The second one
is the search range smaller than a threshold. This might happen when polynomial
has conjugate complex roots. For example, assume polynomial is x^4  + 1. The
sign variation different in interval 1 - alpha to 1 + alpha is always 4. No
matter how small alpha is. However, there should have no real roots. Therefore,
when search range is small enough, we can safely say there is no real root in
this range.

## Page 8

This upper bound here we use Lagrange's bound. Which is the maximum of 1 and sum
of absolute value of ai / an for i from 0 to n-1.

Question here?

## Page 9

Let's go to the continued fraction method. Let's first look at Mobius
transformation. We node it by M(x), which represent ax+b divided by cx+d and ac
- bd not equals to 0.

With mobius transformation, we can map x from 0 to infinity to intervals
represented by a/c and b/d. We can replace x with M(x) in original polynomial to
get a new polynomial P(M(x)). Then the number of  positive roots of A(x) equals
to the number of real root in range a/c and b/d.

## Page 10



