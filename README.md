# About this repository
Deferred correction is a well-established method for incrementally
increasing the order of accuracy of a numerical solution to a set of
ordinary differential equations. Because implementations of deferred
corrections can be pipelined, multi-core computing has increased the
importance of deferred correction methods in practice, especially in
the context of solving initial-value problems. This repository
contains lightweight, flexible MATLAB software to explore variants of
deferred correction methods described in the review article: "Deferred
Correction Methods for Ordinary Differential Equations"
[preprint](http://mathgeek.us/research/papers/dcrev.pdf)

To generate:

- figures 4.1, 4.2 4.3, 4.6, 4.7, run the script: zadunaisky_test.m
- figures 4.4, 4.5, 4.8, 4.9, run the script: idc_test_vector.m and idc_test_vector_work.m
- figure 4.11, run the script: cdc_interpolant/m

