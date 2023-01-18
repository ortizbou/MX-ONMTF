# MX-ONMTF
This is the code for the multiplex community detection method proposed in “Community detection in multiplex networks based on orthogonal nonnegative matrix tri-factorization”\
Authors: Meiby Ortiz-Bouza and Selin Aviyente\
Department of Electrical and Computer Engineering, Michigan State University, MI\
https://arxiv.org/abs/2205.00626v2

## Introduction
In this work, we introduce a new multiplex community detection method that identifies communities that are common across layers as well as those that are unique to each layer. The proposed method, Multiplex Orthogonal Nonnegative Matrix Tri-Factorization, represents the adjacency matrix of each layer as the sum of two low-rank matrix factorizations  corresponding to the common and private communities, respectively. Unlike most of the existing methods which require the number of communities to be pre-determined, the proposed method also introduces a two stage method to determine the number of common and private communities.

<img width="488" alt="Framework_new" src="https://user-images.githubusercontent.com/92049169/212778011-43e63d5d-f825-485e-85d0-78f23a2c0f9e.png">

## Description

Use ```findingk.m``` to find the number of common communities and the number of private and total communities per layer. If the number of communities is known, this step can be skipped.





