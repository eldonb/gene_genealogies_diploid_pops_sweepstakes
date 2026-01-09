\pdfoutput=1
\documentclass[a4paper,10pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
%\usepackage[lf]{Baskervaldx}
%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage{bbm}
%%\usepackage[round,numbers,super]{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\setstretch{1}
\newcommand{\uN}{\ensuremath{\zeta(N)}}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}
\newcommand{\mengi}[1]{\ensuremath{ \left\{#1 \right \} } }
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\IN}{\ensuremath{\mathds{N}} }
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\aths}[1]{\textcolor{violet}{\small \sf #1 }}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  \textsf{\textbf{\@@title}}
\end{flushleft}
\begin{center}
  \textsc{\@@author}
\end{center}
\egroup
} \makeatother
\title{Gene genealogies in diploid 
populations evolving according to sweepstakes reproduction \\ --- approximating $\EE{R_i(n)}$ for the 
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent}
\author{Bjarki Eldon\footnote{ \href{beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{Supported by Deutsche
Forschungsgemeinschaft (DFG) - Projektnummer 273887127
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17
to Wolfgang Stephan; acknowledge funding by the Icelandic Centre of
Research (Rann\'is) through an Icelandic Research Fund
(Ranns\'oknasj\'o{\dh}ur) Grant of Excellence no.\ 185151-051 to Einar
\'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge, Wolfgang
Stephan, and BE;  Start-up module grants through
SPP 1819 with Jere Koskela and Maite
Wilke-Berenguer, and with Iulia Dahmer. \\
\today}\orcidlink{0000-0001-9354-2391} }

\begin{document}

\maketitle
%%\rule{\textwidth}{.5pt}
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


\begin{abstract}
Let $\set{\xi^{n}} \equiv   \set{\xi^{n}(t); t \ge 0}$ denote the
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent for
$0 < \alpha < 2$ and $0 < \gamma \le 1$. Write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \# \set{\xi \in \xi^{n}(t) : \#\xi
= i }dt $, $L(n) \equiv \int_{0}^{\tau(n)} \#\xi^{n}(t)dt $,
$\tau(n) \equiv \inf \set{ t \ge 0 : \#\xi^{n}(t) = 1}$, and
$R_{i}(n) \equiv L_{i}(n)/L(n) $ for $i = 1,2, \ldots, n-1$.  With
this C++ code we estimate the functionals  $\EE{R_{i}(n)}$
for $i = 1,2, \ldots, 2n-1$ when the gene genealogy of $2n$ sampled
gene copies is described by the
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent for
$0 < \alpha < 2$ and $0 < \gamma \le 1$.  The stated coalescent can be
obtained from a population genetics model of a diploid panmictic
population of constant size evolving absent selfing and according to
sweepstakes reproduction, and extends the
$\Omega$-Beta$(2-\alpha,\alpha)$-coalescent \cite{BLS15}.  The
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent is a
continuos-time simultaneous multiple-merger coalescent with an atom at
zero and where the ancestral lineages (blocks of a partition of
$\{1,2,\ldots, 2n\}$) merging in a single group are independently and
uniformly at random and with replacement assigned a label from
$\{1,2,3,4\}$ and the blocks assigned the same label are merged.  Only
one such group of blocks can appear each time, corresponding to the
appearance of one large family involving four parent chromosomes.
\end{abstract}

\tableofcontents


@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline


This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \cite{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \cite{kernighan1988c}
      file.


Compiles on Linux debian 6.12.6-amd64  with ctangle 4.11 and  g++ 14.2 and  GSL 2.8



One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library. 
Using a Makefile can be helpful, naming this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        c++ -Wall -Wextra -pedantic -O3 -march=native -m64 -x c++  iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}



Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v ---leak-check=full ---show-leak-kinds=all  ---leak-resolution=high}

{\tt ---num-callers=40 ---vgdb=full <program start>}

Use {\tt cppcheck} to check the code

{\tt  cppcheck ---enable=all ---language=c++  iguana.c}



To generate estimates on a computer with several CPUs it may be
convenient to put  in a text file ({\tt simfile}):


{\tt ./a.out \$(shuf -i 484433-83230401 -n1) > resout<i>}


for  $i = 1,\ldots, y$ and use  {\tt
parallel}\cite{tange11:_gnu_paral}

{\tt parallel ---gnu -jy :::: ./simfile}



@* {\bf intro}.
\label{sec:intro}

Write $(a)_{n} := a(a-1)\cdots (a - n + 1)$ with $(a)_{0} := 1$, and
for any given event/condition $E$ let $\one{E} := 1$ whenever $E$
occurs/holds, and take $\one{E} = 0$ otherwise; write
$[n] := \mengi{1,2,\ldots,n}$ for any
$n\in \IN := \mengi{1,2,\ldots}$.  Consider a
$\Omega$-$\delta_{0}$-Beta$(\gamma, 2-\alpha,\alpha)$-coalescent on
the partitions of $[2n]$ for $n \in \mengi{2,3,\ldots}$; we interpret
$n$ as the number of diploid individuals and so $2n$ gene copies
sampled.  Write $s = m - k_{1} - \ldots - k_{4}$ for
$k_{1},\ldots, k_{4} \in \{0,2,3, \ldots, m\}$ and
$2 \le k_{1} + \cdots + k_{4} \le m$ for any
$m \in \mengi{2,3, \ldots, n}$.  The
$\Omega$-$\delta_{0}$-Beta$(\gamma, 2-\alpha,\alpha)$-coalescent has
transition rates
\begin{equation}
\label{eq:lambdankstotal}
\lambda_{m;k_{1},\ldots, k_{4};s} =  \binom{m }{k_{1} \ldots k_{4} \quad  s} \frac{1}{\prod_{j=2}^{m} \svigi{\sum_{i} \one{k_{i} = j }}!   } \lambda_{m,k_{1},\ldots,k_{4},s}^{\prime}
\end{equation}
where  $\lambda_{m;k_{1},\ldots, k_{4};s}$ is the total rate of merging blocks in groups of sizes $k_{1},\ldots, k_{4}$ and $r = \sum_{j}\one{k_{j} \ge 2}$ 
\begin{equation}
\label{eq:lambdanks}
\lambda_{m,k_{1},\ldots,k_{4},s}^{\prime} =  \one{\sum_{j}k_{j}=2}\frac{C_{\kappa}}{C_{\alpha,\gamma}} + \frac{\alpha c}{C_{\alpha,\gamma} \mathbbm m ^{\alpha} }\sum_{\ell=0}^{s\wedge (4-r)} \binom{s}{\ell} \frac{(4)_{r+\ell} }{4^{k+\ell}} B(\gamma, k +\ell  - \alpha, m-k -\ell +\alpha)
\end{equation}
where
$B(p,a,b) = \int_{0}^{1} \one{ 0 < u \le p }u^{a-1}(1-u)^{b-1}du$ for
$0 < p \le 1$ and $a,b > 0$ and
\begin{subequations}
\begin{align}
\label{eq:3}
 \mathbbm {m}  & =   (2 + (1 + 2^{1-\kappa})/(\kappa - 1))/2  \\
\label{eq:gammas}
\gamma & =   \one{ \frac{\uN}{N} \to K } \frac{ K }{K + \mathbbm m } + \one{ \frac{\uN}{N} \to \infty } \\
\label{eq:2}
C_{\alpha,\gamma} & =  \frac 14 C_{\kappa} +   \frac{ \alpha c }{ 4\mathbbm{m}^{\alpha} }B(\gamma, 2-\alpha, \alpha) \\
\label{eq:1}
C_{\kappa} & =  \one{\kappa = 2} \frac{2}{\mathbbm  m^{2}} +   \one{\kappa > 2}  \frac{2}{\mathbbm m^{2}} \frac{c_{\kappa} }{2^{\kappa}(\kappa - 2)(\kappa - 1) }
\end{align}
\end{subequations}
where in \eqref{eq:1} $\kappa + 2 < c_{\kappa} < \kappa^{2}$ when
$\kappa > 2$.  In \eqref{eq:lambdanks}  $c > 0$ and  $0 < \alpha < 2$. 


In \eqref{eq:lambdanks} we take $0 < \alpha < 2 $ with the
understanding that the population model when $0 < \alpha < 1$ is
different from the one when $1 \le \alpha < 2$.  Let $X$ be a positive
integer-valued random variable with law
\begin{equation}
\label{eq:lawX}
\prb{X = k}      =   C\svigi{k^{-\alpha}  -  (1+k)^{-\alpha}}
\end{equation}
for $k \in \mengi{2,3,\ldots, \uN}$ where $C$ is such that
$\prb{2 \le X \le \uN} = 1$; $X$ is the random number of potential
offspring of an arbitrary parent pair.  In any given generation it is
assumed that the $2N$ current individuals randomly form parent pairs,
and the pairs then independently produce potential offspring according
to \eqref{eq:lawX}.  Out of the at least $2N$ potential offspring so
generated we sample $2N$ uniformly at random and without replacement
to survive and replace the current individuals.  Write
$X \vartriangleright \mathds L(a,\uN)$ when $X$ is given law
\eqref{eq:lawX} with $a$ and $\uN$ as given each time.  Write
$X_{1}(g), \ldots, X_{N}(g)$ for the random number of potential
offspring produced in generation $g$, and $\svigi{U_{g}}_{g}$ for a
sequence of i.i.d.\ random uniforms and $(\varepsilon_{N})_{N}$ where $0 < \varepsilon_{N} < 1$.
 Suppose  $\kappa,  \uN \ge 2$    fixed and  $X_{1}(g), \ldots, X_{N}(g)$ are iid copies of $X(g)$ where 
\begin{equation}
\label{eq:randall}
X(g) \vartriangleright \mathds L\svigi{ \one{U_{g} \le \varepsilon_{N}} \alpha + \one{ U_{g} > \varepsilon_{N}} \kappa, \uN }
\end{equation}
when $1 \le \alpha < 2$; when $0 < \alpha < 1$ it holds that 
\begin{equation}
\label{eq:randone}
X_{i}(g) \vartriangleright \mathds L\svigi{ \one{U_{g} \le \varepsilon_{N}} \alpha + \one{ U_{g} > \varepsilon_{N}} \kappa, \uN }, \quad X_{j\neq i}(g)  \vartriangleright \mathds L\svigi{  \kappa, \uN}
\end{equation}
In \eqref{eq:randone} the index $i$ is picked uniformly at random from
$[N]$, and $X_{j}(g) \vartriangleright L\svigi{ \kappa, \uN} $ for all
$j \neq i$.  Both \eqref{eq:randall} and \eqref{eq:randone} lead to
\eqref{eq:lambdanks}, so that when we use \eqref{eq:lambdanks} with
$0<\alpha < 2$ it is with the understanding that \eqref{eq:randone}
holds when $0<\alpha < 1$, and when $1 \le \alpha < 2$
\eqref{eq:randall} is in force.  In \eqref{eq:randall} the
$X_{1},\ldots, X_{N}$ are i.i.d.; they are independent but may not
always be identically distributed when \eqref{eq:randone} holds.


Let $\mengi{\xi^{n}} \equiv \mengi{\xi^{n}(t); t \ge 0}$ be the coalescent
on the partitions of $[2n]$ with transitions rates
\eqref{eq:lambdanks} and   $\xi^{n}(0) =  \mengi{ \mengi{1}, \ldots, \mengi{2n}   }$.    Define the functionals, where  $\#A$ is
the cardinality of a given  set $A$, 
\begin{equation}
\label{eq:L}
L_{i}(n) :=  \int_{0}^{\tau(n)}  \#\mengi{\xi \in \xi^{n}(t) : \#\xi = i}dt, \quad  L(n) :=   \int_{0}^{\tau(n)}  \#\xi^{n}(t)dt
\end{equation}
where $i\in [2n-1]$ and
$\tau(n) := \inf\mengi{ t \ge 0 : \# \xi^{n}(t) = 1 }$.  Then
$L_{i}(n)$ are the random length of branches supporting $i$ leaves
(gene copies), and $L(n) = L_{1}(n) + \cdots + L_{2n-1}(n)$ is the
random total tree length\cite{BLS15}.  Define, for
$n\in \{2,3,\ldots\}$,
\begin{equation}
\label{eq:R}
R_{i}(n) :=   \frac{L_{i}(n) }{L(n) }, \quad i = 1,2,\ldots, 2n-1
\end{equation}
and $R_{i}(n)$  is well defined since $L(n) > 0$ almost surely.  We will
estimate $\EE{R_{i}(n)}$ when $\mengi{\xi^{n}}$ has transition rates
\eqref{eq:lambdankstotal} with $\lambda_{n;k_{1}, \ldots, k_{r};s}$ as
in \eqref{eq:lambdanks}.


In \S~\ref{sec:code} we briefly summarize the algorithm; the code
follows in \S~\ref{sec:includes}--\S~\ref{sec:main}; we conclude in
\S~\ref{sec:concl}.



@* {\bf code}.
\label{sec:code}

 Let 
\begin{displaymath}
\mathds K_{m} := \mengi{ (k_{1}, \cdots, k_{4}) : k_{j} \in \{0,2,3,\ldots, m\}, \quad    k_1 \ge k_2 \ge k_{3}  \ge k_4, \quad  2 \le   k_{1} +
\cdots + k_{4} \le m}
\end{displaymath}
be the set of all possible ordered merger sizes when there are $m$
blocks and $j_{m} : \mathds K_{m} \to [\# \mathds K_{m}]$ be a
bijection assigning indexes to the mergers ordered such that (see \S~\ref{sec:lambdan} for the ordering)
\begin{displaymath}
1 = j_{m}( (2,0,0,0))  < j_{m}((2,2,0,0)) < \ldots < j_{m}((m,0,0,0)) =  \#\mathds K_{m}.
\end{displaymath}
  Let $s = m - \sum_{j}k_{j}$ and   $\lambda_{m} := \sum_{(k_{1},\ldots, k_{4}) \in \mathds K_{m} } \lambda_{m;k_{1}, \ldots, k_{4};s}$ denote
the total jump rate out of $m$ blocks (recall \eqref{eq:lambdankstotal}).    Let 
$F_{m}: [\#\mathds K_{m}] \to [0,1]$ denote the cumulative mass function 
\begin{displaymath}
  F_{m}(\ell) = \frac{1}{\lambda_{m}} \sum_{i=1}^{\ell} { \lambda_{m;j_{m}^{-1}(i);s}  } 
\end{displaymath}
We use $F_m$ to sample merger sizes. Let $U$ denote a standard random
uniform. Then, given a sample of $U$, $j_{m}^{-1}(j^{*}) \in \mathds
K_{m}$ where $j^{*} := \inf \mengi{j \in [\# \mathds K_{m} ] : U \le
F_{m}(j) }$.
 
Let $\svigi{\ell_{i}(n), \ldots, \ell_{n-1}(n) }$ denote the realised
branch lengths \eqref{eq:L}, and $\svigi{b_{1},\ldots, b_{m}}$ the
current block sizes where $b_{j}\in [2n]$ and $\sum_{j}b_{j} = 2n$
(the sample consists of  $n$ diploid individuals and so $2n$ gene copies). 

\begin{enumerate}
\item $\svigi{ \mathfrak  r_{1}(n), \ldots, \mathfrak  r_{2n-1}(n)} \leftarrow (0,\ldots, 0)$
\item compute $\lambda_{m}$ for $m = 2,3,\ldots, 2n$ given parameter values in \S~\ref{sec:constants};   \S~\ref{sec:lambdanalln} 
\item for each of $M$ experiments   \S~\ref{sec:estimate}
\begin{enumerate}
\item $\svigi{\ell_{1}(n), \ldots, \ell_{2n-1}(n)} \leftarrow (0,\ldots, 0) $
\item $m \leftarrow 2n$ where $m$ is current number of blocks and $2n$ is sample size 
\item $\svigi{b_{1}, \ldots, b_{2n} } \leftarrow \svigi{1,\ldots, 1} $ 
\item {\bf while} $m > 1$: \S~\ref{sec:oneexperiment} 
\begin{enumerate}
\item sample a random exponential $t$ with rate $\lambda_{m}$ 
\item  $\ell_{b}(n) \leftarrow \ell_{b}(n) + t $ for $b =  b_{1}, \ldots, b_{m}$ \S~\ref{sec:updatebranchlengths}
\item sample merger sizes  $j_{m}^{-1}(j^{*}) \in \mathds K_{m}$ where   $j^{*} :=  \inf \mengi{j \in [\# \mathds K_{m} ] : U \le F_{m}(j) }$ \S~\ref{sec:samplemergersizes} 
\item merge blocks according to sampled merger sizes $j_{m}^{-1}(j^{*})  \in \mathds K_{m}$  \S~\ref{sec:mergeblocks}
\item $m\leftarrow m - \sum_{i}k_{i} + \sum_{j}\one{k_{j} \ge 2}$ where $(k_{1}, \ldots, k_{4}) = j_{m}^{-1}(j^{*})$
\end{enumerate}
\item  $\mathfrak r_{i}(n) \leftarrow  \mathfrak r_{i}(n) +  \ell_{i}(n)/\sum_{j}\ell_{j}(n)$  for $i \in [2n-1]$ \S~\ref{sec:updaterelativebranchlengths} 
\end{enumerate}
\item return  an estimate $\overline \varrho_{i}(n) =   (1/M)\mathfrak r_{i}(n) $   of $\EE{R_{i}(n)}$ 
\end{enumerate}



@*1 {\bf includes}.  
\label{sec:includes}


the included libraries; we use the GSL library


@<includes@>=@#
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <chrono>
#include <limits>
#include <cfloat>
#include <assert.h>
#include <math.h>
#include <fenv.h>
#include <errno.h>
#include <unistd.h>
#include <omp.h>
#include <sys/param.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include "xindbetacoal.h"



@*1 {\bf constants}.
\label{sec:constants}

the parameter values 

@<constants@>=@#
const double  dblepsilon = 2.2204460492503131e-16 ;
/* \newline $0 < \alpha < 2$ */
const double CONST_ALPHA = 0.5 ;
/* \newline $\gamma$ \eqref{eq:gammas} */
const double CONST_GAMMA = 0.5 ;
const double CONST_C     = 1.  ;
/* \newline $\kappa \ge 2$ */
const double c_kappa = 2.  ;
/* \newline $\mathbbm m =   (2 +  (1 + 2^{1-\kappa})/(\kappa - 1))/2  $ \eqref{eq:3} */
const double CONST_Minf  = (2. +  (1 + pow(2., 1. - c_kappa))/( c_kappa - 1.))/2. ;
/* \newline $(\kappa + 2 + \kappa^{2})/(2^{1 + \kappa }(\kappa - 2)(\kappa - 1)) =   c_{\kappa}$ */
const double ck = (2. + c_kappa + pow(c_kappa, 2.))/( pow(2.,1. + c_kappa) * (c_kappa-2.)*(c_kappa - 1.) ) ;
/* \newline $C_{\kappa} =  2(\one{\kappa = 2} + \one{\kappa > 2}c_{\kappa}) /\mathbbm m^{2}$ \eqref{eq:1}  */
const double CONST_Ckappa =  2.*( fmax( c_kappa - 2., std::numeric_limits<double>::epsilon() ) > DBL_EPSILON ? ck : 1.) / pow( CONST_Minf, 2.);
/* \newline $C_{\alpha,\gamma}$ \eqref{eq:2} */
const double CONST_Calphagamma =  ((CONST_ALPHA * CONST_C * gsl_sf_beta(2-CONST_ALPHA, CONST_ALPHA) * ( fmax(1. - CONST_GAMMA, DBL_EPSILON) > DBL_EPSILON ? gsl_sf_beta_inc(2. - CONST_ALPHA, CONST_ALPHA, CONST_GAMMA) : 1.)/pow(CONST_Minf,CONST_ALPHA)) + CONST_Ckappa)/4. ;
/* \newline sample size */
const std::size_t CONST_SAMPLE_SIZE = 30 ;
/* \newline number of experiments */
const int CONST_EXPERIMENTS = 25e4; 



@*1 {\bf random number generators}. 
\label{sec:rngs}

define the STL |rng()| and GSL |rngtype|  random number generators 

@<rngs@>=@#
/* \newline \aths{ random seed generator for the STL random number generator} */
  std::random_device randomseed;
  /* \newline \aths{  Standard mersenne twister  random number engine seeded with |randomseed()|}  */
  std::mt19937_64 rng(randomseed());


/* \newline  \aths{ define the  GSL random number generator |rngtype|} */
gsl_rng * rngtype;

static  void setup_rng(  const unsigned long s)
  {
    const gsl_rng_type *T ; 
    gsl_rng_env_setup(); 
    T = gsl_rng_default ;
    rngtype = gsl_rng_alloc(T);
    gsl_rng_set( rngtype,  s) ;
  }


@*1 {\bf the exponential function}.
\label{sec:expfunction}

compute $e^{x}$ for the given $x$ with error checking 

@<exponentialfunction@>=@#
static long double veldi( const long double & v )
{
  feclearexcept(FE_ALL_EXCEPT);

  const long double svar =  expl( v ) ;

  assert( fetestexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW) == 0 );

  return ( svar );
}




@*1 {\bf the beta function}.
\label{sec:betafunc}

compute the logarithm of the  (incomplete) beta function
$\int_{0}^{1}\one{0 < u \le x}u^{a-1}(1-u)^{b-1} du$ using the Gauss
hypergeometric function
\begin{displaymath}
B(x,a,b) = \frac 1a x^{a}(1-x)^{b}F(a+b, 1, a+1, x)
\end{displaymath}

@<betafunction@>=@#
static long double betafunc( const double& a, const  double&  b )
{
  /* \newline \aths{  the GSL incomplete beta function is normalised by the complete beta function }  */
  
  assert( a > DBL_EPSILON);
  assert( b > DBL_EPSILON);
  
  /* \newline  \aths{ the standard way would be | gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) | } */
  /* \newline \aths{ return  the logarithm  of the beta function  as $\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b)$ } */
  const long  double f = static_cast<long double>( (  1. - CONST_GAMMA > DBL_EPSILON ? gsl_sf_hyperg_2F1(a + b, 1., a+1., CONST_GAMMA) : 1.) );
  
  assert( f > DBL_EPSILON );
  /* \newline \aths{  return $\one{x < 1}( \log f  +a\log x + b\log(1-x) - \log a) + \one{x = 1} (\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b)) $} */
  return( 1. - CONST_GAMMA  >  DBL_EPSILON ?  (logl( f ) +  ( static_cast<long double>( (a*log(CONST_GAMMA)) + (b * log(1. - CONST_GAMMA)) - log(a) ) ) ) : (lgammal( static_cast<long double>(a)) + lgammal( static_cast<long double>(b) ) - lgammal( static_cast<long double>(a + b) ) ) );
}



@*1 {\bf the counts constant}.
\label{sec:binomconst}

given merger sizes $k_{1}, \ldots, k_{r}$ compute
$\prod_{j=2}^{m} \svigi{\sum_{i} \one{k_{i} = j} }!  $  \eqref{eq:lambdankstotal}

@<countconstant@>=@#
static long double numbercollisions( const std::vector<std::size_t>& __k )
{
 
  assert( __k[0] > 1);
  /* \newline get the number of (simultaneous) mergers  $|r|  \in [4]$ */
  const int  r = (__k[3]> 1 ? 4 : (__k[2] > 1 ? 3 : (__k[1] > 1 ? 2 : 1))) ;

  int  l {} ;
  switch( r ){
  case 1 : {
    l = 1;  break ; }
  case 2 : {
    assert( __k[1] <= __k[0]); 
    l = (__k[0] == __k[1] ? 2 : 1);
    break ; }
  case 3 : {
    assert( __k[2] <= __k[1]);
    assert( __k[1] <= __k[0]);
    l = ( __k[0] == __k[2] ? 6 : ( __k[0] == __k[1] ? 2 : ( __k[1] == __k[2] ? 2 : 1)));
    break ; }
  case 4 : {
    assert( __k[3] <= __k[2]);
    assert( __k[2] <= __k[1]);
    assert( __k[1] <= __k[0]);
    l = ( __k[0] ==  __k[3] ? 24 : ( __k[0] == __k[2] ? 6 :  ( ( __k[0] == __k[1] ? ( __k[2] == __k[3] ? 4 : 2) : ( __k[1] == __k[3] ? 6 : ( __k[1] == __k[2] ? 2 : ( __k[2] == __k[3] ? 2 :1)))))));
    break ; }
  default : break ; }


 assert( l > 0);
 /* \newline  return   $\prod_{j=1}^{m} \svigi{\sum_{i} \one{k_{i} = j} }!  $ */
  return static_cast<long double>( l)  ;
}


@*1 {\bf the falling factorial}.
\label{sec:fallingfact}

compute the logarithm of the falling factorial
$(4)_{m} := \prod_{j=0}^{m-1} (4-j)$ with $(4)_{0} := 1$

@<falling@>=@#
static long double ff( const std::size_t &  m )
{

 assert( m <  5 ); 
  return logl( static_cast<long double>(m > 0 ? (m < 2 ? 4 : ( m < 3 ? 12 : 24) ) : 1) ) ;

}


@*1 {\bf the multinomial constant $\tbinom{m }{k_{1}\ldots k_{r}\quad s}$ }. 
\label{sec:multinomialconst}


return  the logarithm of the   multinomial constant  $\tbinom{m }{k_{1}\ldots k_{r}\quad s}$ in \eqref{eq:lambdankstotal} as $\log \Gamma (1+m) - \log\Gamma(1+k_{1}) - \one{k_{2} \ge  2}\log\Gamma(1+k_{2}) -  \one{k_{3} \ge  2}\log\Gamma(1+k_{3}) -    \one{k_{4} \ge  2}\log\Gamma(1+k_{4}) - \log\Gamma(1+m-\sum_{i}k_{i}) $

@<multinomial@>=@#
static long double multinomialconstant( const std::size_t& m,   const std::vector<std::size_t>& v__k )
{
assert( m > 1);
     assert( v__k[0] > 1);

 long double svar = lgammal( static_cast<long double>(m+1)) ;
  svar -=  lgammal(  static_cast<long double> ( v__k[1] + 1) );
  svar -=  ( v__k[2] > 1 ? lgammal(  static_cast<long double> ( v__k[2] + 1 )) : (long double)0.) ;
  svar -=  ( v__k[3] > 1 ? lgammal(  static_cast<long double> ( v__k[3] + 1 )) : (long double)0.) ;
  svar -=  ( v__k[4] > 1 ? lgammal(  static_cast<long double> ( v__k[4] + 1 )) : (long double)0.) ;
  svar -=  lgammal(  static_cast<long double> (1+ m -  v__k[1] - v__k[2] - v__k[3] - v__k[4]) )  ;

return svar ;
}





@*1 {\bf the total merger rate $\lambda_{m;k_{1},\ldots, k_{r};s}$ \eqref{eq:lambdankstotal}}.
\label{sec:lambdankstotal}


given merger sizes the total merger rate
$\lambda_{m;k_{1},\ldots, k_{r};s}$ \eqref{eq:lambdankstotal}

@< $\lambda_{m;k_{1},\ldots, k_{r};s}$ \eqref{eq:lambdankstotal} @>=@#
static long double  lambdankstotal( const std::size_t & m,    const std::vector<std::size_t>&   __k )
{

  @#

   assert( __k[0] > 1);
 
  assert( __k[3] <= __k[2] );
  assert( __k[2] <= __k[1] );
  assert( __k[1] <= __k[0] );
  const std::size_t r = 1 + (__k[1] > 1 ? 1 : 0) +  (__k[2] > 1 ? 1 : 0) +  (__k[3] > 1 ? 1 : 0) ;
  assert(r > 0);
  const std::size_t  k =   __k[0] + __k[1] + __k[2] + __k[3]  ;
  assert( k > 1);

  long double x {} ;
  long double y {} ;
  for( std::size_t ell = 0 ; ell <= MIN(m-k, 4-r); ++ell){
  /* \newline  \aths{ $\log \binom{m}{k_{1}\ldots k_{4}\quad s }$  \S~\ref{sec:multinomialconst} } */
    y =  multinomialconstant(m, __k) ;
    /* \newline \aths{  $\log \binom{s}{\ell}$ from  \eqref{eq:lambdanks}}  */
    y +=  static_cast<long double>(gsl_sf_lnchoose(m-k, ell)) ;
    /* \newline $(4)_{r+\ell}$ from   \eqref{eq:lambdanks}; \S~\ref{sec:fallingfact} */
    y +=  ff( r + ell) ;
    /* \newline $B(\gamma, \sum_{j}k_{j} + {\ell} - \alpha, m-\sum_{j}k_{j}-\ell + \alpha)$ from   \eqref{eq:lambdanks}; \S~\ref{sec:betafunc} */
     y +=  betafunc(  static_cast<double>(k + ell) - CONST_ALPHA, static_cast<double>(m-k-ell) + CONST_ALPHA) ;
     /* \newline $\log 4^{\sum_{j} k_{j} + \ell}$ from    \eqref{eq:lambdanks} */
    y -=  (static_cast<long double>(k+ell) * logl( static_cast<long double>(4) ) ) ;
    /* \newline  $\exp\left( \log \binom{s}{\ell} + \log (4)_{r+\ell} + \log B(\gamma,  \sum_{j}k_{j} + {\ell} - \alpha, m-\sum_{j}k_{j}-\ell + \alpha) - \log 4^{ \sum_{j} k_{j} + \ell} \right)$ recalling   \eqref{eq:lambdanks};   \S~\ref{sec:expfunction} */
    x += veldi(y) ; }

/* \newline   $\prod_{j=1}^{m} \svigi{\sum_{i} \one{k_{i} = j} }! $ from  \eqref{eq:lambdankstotal};  \S~\ref{sec:binomconst} */
  x /= numbercollisions(__k) ;


 /* \newline multiply with  $\alpha c/ \mathbbm m ^{\alpha}$ from  \eqref{eq:lambdanks} */
  x *=  (CONST_ALPHA * CONST_C / ( powl( static_cast<long double>(CONST_Minf), static_cast<long double>(CONST_ALPHA)) ) ) ;

/* \newline add $C_{\kappa} m(m-1)/2$ from  \eqref{eq:lambdanks} when $\sum_{j}k_{j} = 2$   */
  x += (k < 3 ? static_cast<long double>( (m*(m-1))/2 ) * static_cast<long double>( CONST_Ckappa) : static_cast<long double>(0) ) ;
  x /= static_cast<long double>( CONST_Calphagamma ) ;
  
  return ( x ) ;

}


@*1 {\bf the total jump rate out of $m$ blocks}.
\label{sec:lambdan}

the total jump rate $\lambda_{m} = \sum_{(k_{1},\ldots, k_{4})\in
\mathcal K_{m}} \lambda_{m;k_{1},\ldots, k_{4}; s}$ out of $m$ blocks
where $ \lambda_{m;k_{1},\ldots, k_{r}; s}$ as in
\eqref{eq:lambdankstotal}.  Need to sum over all possible ordered
merger sizes $\bold{k} \in \mathcal K_{m} $ given number of blocks $m$


@<lambdan@>=@#
static void lambdan( const std::size_t & m, std::vector<long double>& v__lambdan )
{
  /*  \newline       |m| is current number of blocks */
  assert( m > 1);
  /* \newline $\_\_k$  will contain the merger sizes */ 
  std::vector<std::size_t> __k (4); 
  
   for( std::size_t i = 2 ; i <= m ; ++i){
    __k = {i, 0, 0, 0};
    assert( std::accumulate( __k.begin(), __k.end(), 0) <= static_cast<int>(m) );
    /* \newline \S~\ref{sec:lambdankstotal} */
    v__lambdan[m] +=  lambdankstotal(m, __k)   ;

    if( MIN(i, m-i) > 1 ){
      for( std::size_t j = 2 ; j <= MIN( i, m-i); ++j)
	{
	  __k = {i,j, 0, 0}; 
	  assert( __k[1] <= __k[0]); 
	  assert( std::accumulate( __k.begin(), __k.end(),0) <= static_cast<int>(m) );
	  
	  v__lambdan[m] +=  lambdankstotal(m, __k) ;
	  
	  if( MIN(j, m-i-j) > 1){
	    for( std::size_t k = 2; k <=  MIN(j,m-i-j); ++k){
	      assert( k <= j);
	      __k = {i, j, k, 0} ;

	      assert( std::accumulate( __k.begin(), __k.end(),0) <= static_cast<int>(m) );
	      
	      v__lambdan[m] +=  lambdankstotal(m, __k) ;
	      
	      if( MIN(k, m-i-j-k) > 1){
		for( std::size_t l  = 2; l <= MIN(k, m-i-j-k); ++l){
		  
		  __k = {i, j, k, l} ;

		  
		  assert( std::accumulate( __k.begin(), __k.end(),0) <= static_cast<int>(m) );
		 v__lambdan[m] +=  lambdankstotal(m, __k);
		   ;}}}}}}}
}


@*1 {\bf $\lambda_{m}$ for $m = 2,3,\ldots, n$}. 
\label{sec:lambdanalln}

compute the total jump rate $\lambda_{m}$ for $m = 2,3,..., n $ where
$n$ is sample size; for each $m = 2, 3, \ldots, n$ add up
$\lambda_{m;k_{1}, \ldots, k_{r};s}$ from \eqref{eq:lambdankstotal}
over the ordered merger sizes  
 
@<total jump rate all $m$@>=@#
static void lambdamallm( std::vector<long double>& v__lambdam)
{
 

  /* \newline  |__lambdam| stores the $\lambda_{m}$ values;   |CONST_SAMPLE_SIZE| \S~\ref{sec:constants} */
  for( std::size_t n = 2; n <= CONST_SAMPLE_SIZE; ++n){
  /* \newline |lambdan| \S~\ref{sec:lambdan} */
     lambdan( n,  v__lambdam ); }
}


@*1 {\bf sample merger sizes}. 
\label{sec:samplemergersizes}


sample merger sizes by going through the ordered mergers until
$u \le F(\bold k)$ where $u$ is a random uniform and $F$ is the sum of
the merger probabilities up to and including $\bold k$

@<sample $\bold k$@>=@#
static void samplemergersizes( const std::size_t& m,  const std::vector<long double>& __lambdan, std::vector<std::size_t>& __k )
{
  @#

  assert( __k.size() == 4 );
  std::size_t i = 1 ;
  std::size_t j {} ;
  std::size_t k {} ;
  std::size_t l {} ;

  long double F {} ;
  const double u = gsl_rng_uniform(rngtype);
  while( (i <= m) && (u > F) ){
    ++i ;
    __k = { i, 0, 0, 0} ;
    /* \newline  \aths{ |lambdankstotal| \S~\ref{sec:lambdankstotal}} */
     F += lambdankstotal(m, __k) / __lambdan[m] ;
    
    if( (MIN( i, m-i) > 1) && (u > F) ){
      j = 1 ;
      while( (j <=  MIN( i, m-i)) && (u > F) ){
	++ j;
	__k = {i, j, 0, 0}; 
	  F += lambdankstotal(m, __k) / __lambdan[m] ;
	
	if( (MIN(j, m-i-j) > 1) && (u > F) ){
	  k = 1 ;
	  while( (k <= MIN(j, m-i -j) ) && (u > F) ){
	    ++k ;
	    __k = {i, j, k, 0}; 
	      F += lambdankstotal(m, __k) / __lambdan[m] ;
	    
	    if( (MIN( k, m - i - j - k) > 1) && (u > F)){
	      l = 1 ; 
	      while( (l <= MIN(k,m-i-j-k)) && (u > F)){
		++ l;
		 __k = {i, j, k, l}; 
		   F += lambdankstotal(m, __k) / __lambdan[m];  }}}}}}}
		 
  assert( static_cast<std::size_t>( accumulate( __k.begin(), __k.end(), 0)) <= m) ;
}


@*1 {\bf merge blocks}. 
\label{sec:mergeblocks}

merge blocks and record new block sizes given merger sizes

@<merge blocks given merger sizes@>=@#
static void  mergeblocks( const std::vector<std::size_t>& __k, std::vector<std::size_t>& __blocks )
{

  /* \newline \aths{  |__k|  is merger sizes; |__blocks| the current block sizes} */
  assert( __k[0] > 1 );
    std::shuffle(   __blocks.begin(), __blocks.end(), rng );
        
  const std::size_t r = ( __k[3] > 1 ? 4 : (__k[2] > 1 ? 3 : ( __k[1] > 1 ? 2 : 1))) ;
  std::vector< std::size_t> newblocks (r);
  const std::size_t m = __blocks.size() ;

  for( std::size_t i = 0 ; i < r ; ++i){
    assert(  __k[i] <=  __blocks.size() );
    newblocks[i] = static_cast<std::size_t>( std::accumulate( __blocks.rbegin(), __blocks.rbegin() +  __k[i], 0));

    assert( newblocks[i] > 1);
    __blocks.resize( __blocks.size() - __k[i]);}
  __blocks.insert( __blocks.end(), newblocks.begin(), newblocks.end() );

  assert( __blocks.size() == ( m - (__k[0] + __k[1] + __k[2] + __k[3]) + r ) ) ;
}


@*1 {\bf  update branch lengths $\ell_{i}(n)$}. 
\label{sec:updatebranchlengths}


given current block sizes $\_\_b$ update branch lengths
$\ell_{b}(n) \leftarrow t + \ell_{b}(n)$ for $b = b_{1},\ldots, b_{m}$
where $t$ is a random exponential with rate $\lambda_{m}$ 

@<update $\ell_{i}(n)$ @>=@#
static void updatebranchlengths( const double &t, const std::vector<std::size_t>& __b,  std::vector<double>& __l  )
{

@#

  for (const auto &b: __b){
    assert( b > 0);
    assert( b < CONST_SAMPLE_SIZE);
    __l[0] += t ;
    __l[b] += t ; }
}


@*1 {\bf update relative branch lengths $r_{i}(n)$}.
\label{sec:updaterelativebranchlengths}

given a realisation of branch lengths update estimate
$\overline \varrho _{i}(n)$ of mean relative branch lengths

@<update $\overline \varrho_{i}(n)$@>=@#
static void updaterelativebranchlengths( const std::vector<double>& __l, std::vector<double>& __r)
{
 
 /* \newline \aths{ $\_\_l$ is vector of branch  lengths; $\_\_r$ is vector of estimates  $\overline \varrho_{i}(n) $ } */
  assert( __l[0] > 0);

  const double d = __l[0] ;
  std::transform( __l.begin(), __l.end(), __r.begin(), __r.begin(), [&d]( const auto &x, const auto &y){ return y +  (x/d);});    
}

@*1 {\bf one experiment}.
\label{sec:oneexperiment}

generate one realisation of branch lengths - one experiment

@<one experiment@>=@#
static void one( std::vector<double>&  __errs, const std::vector<long double>& __lambdam )
{

 @#

  std::vector<double> __ells (CONST_SAMPLE_SIZE) ;

 @#

  std::vector<std::size_t> __bees (CONST_SAMPLE_SIZE, 1) ;

 @#

  std::vector<std::size_t> __mergersizes (4); 

  double t {} ;
  while( __bees.size() > 1){

/* \newline \aths{ sample a random exponential with rate $\lambda_{m}$ } */
    t = gsl_ran_exponential( rngtype, 1. / __lambdam[ __bees.size() ] ) ;

/* \newline \S~\ref{sec:updatebranchlengths} */
    updatebranchlengths( t,  __bees,  __ells);

/* \newline \aths{ sample merger sizes \S~\ref{sec:samplemergersizes} } */
    samplemergersizes( __bees.size(), __lambdam, __mergersizes);

/* \newline \aths{ given merger sizes merge blocks \S~\ref{sec:mergeblocks}} */
    mergeblocks( __mergersizes,  __bees); }

  /* \newline  \aths{ given branch lengths update relative branch lengths \S~\ref{sec:updaterelativebranchlengths}} */
  updaterelativebranchlengths( __ells, __errs); 
}



@*1 {\bf estimate $\EE{R_{i}(n)}$}.
\label{sec:estimate}


obtain an estimate $\overline \varrho_{i}(n)$ of $\EE{R_{i}(n)}$
\eqref{eq:R} for the given number of experiments defined in
\S~\ref{sec:constants}

@<estimate@>=@#
static void estimate()
{
/* \newline \S~\ref{sec:constants} */
  int M = CONST_EXPERIMENTS + 1 ;

  std::vector<long double> __lambdam (1 + CONST_SAMPLE_SIZE, 0 );
  std::vector<double> __varrho (CONST_SAMPLE_SIZE, 0);
  /* \newline  \aths{  get the total rate \S~\ref{sec:lambdanalln} } */
  lambdamallm(  __lambdam ) ;

  /* \newline \aths{ |one|  \S~\ref{sec:oneexperiment} } */
  while( --M > 0){ one( __varrho, __lambdam);}

  for( const auto &v : __varrho ){ std::cout << v << '\n' ;}
}





@*1 {\bf the main module}.
\label{sec:main}

the |main| function

@C 

/* \newline \S~\ref{sec:includes} */
@<includes@>@#
/* \newline \S~\ref{sec:constants} */
@<constants@>@#
/* \newline \S~\ref{sec:rngs} */
@<rngs@>@#
/* \newline \S~\ref{sec:expfunction} */
@<exponentialfunction@>@#
/* \newline \S~\ref{sec:betafunc} */
@<betafunction@>@#
/* \newline \S~\ref{sec:binomconst} */
@<countconstant@>@#
/* \newline \S~\ref{sec:fallingfact} */
@<falling@>@#
/* \newline \S~\ref{sec:multinomialconst} */
@<multinomial@>@#
/* \newline \S~\ref{sec:lambdankstotal} */
@< $\lambda_{m;k_{1},\ldots, k_{r};s}$ \eqref{eq:lambdankstotal} @>@#
/* \newline \S~\ref{sec:lambdan} */
@<lambdan@>@#
/* \newline \S~\ref{sec:lambdanalln} */
@<total jump rate all $m$@>@#
/* \newline \S~\ref{sec:samplemergersizes} */
@<sample $\bold k$@>@#
/* \newline \S~\ref{sec:mergeblocks} */
@<merge blocks given merger sizes@>@#
/* \newline \S~\ref{sec:updatebranchlengths} */
@<update $\ell_{i}(n)$ @>@#
/* \newline \S~\ref{sec:updaterelativebranchlengths} */
@<update $\overline \varrho_{i}(n)$@>@#
/* \newline \S~\ref{sec:oneexperiment} */
@<one experiment@>@#
/* \newline \S~\ref{sec:estimate} */
@<estimate@>@#

int main(int argc, char * argv[])
{
  
 /*  \newline  \aths{ |setup_rng| \S~\ref{sec:rngs}} */
 setup_rng( static_cast< std::size_t>( atoi(argv[1])) ) ;

 /* \newline \aths{ |estimate| \S~\ref{sec:estimate} } */
  estimate() ;

return GSL_SUCCESS ;
}



@* {\bf conclusions and bibliography}. 
\label{sec:concl}


We estimate mean relative branch lengths $\EE{R_{i}(n)}$ as predicted
by the
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent. The
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent can be
derived from a model of a diploid panmictic population of constant
size evolving absent selfing, in a random environment, and according
to randomly increased recruitment.  The
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent
extends the $\Omega$-Beta$(2-\alpha,\alpha)$-coalescent derived from
diploid panmictic populations \cite{BLS15} and based on
\cite{schweinsberg03}. The
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent
extends the $\Omega$-Beta$(2-\alpha,\alpha)$-coalescent by {\it (i)}
including an atom at zero, {\it (ii)} extending the range of $\alpha$
from $[1,2)$ to $(0,2)$, and {\it (iii)} including a truncation
$0 < \gamma \le 1$.  Moreover {\it (iv)}, the unit of time of the
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent is
proportional to (at least) $N/\log N$ generations (where $2N$ is the
population size); the unit of time for the
$\Omega$-Beta$(2-\alpha,\alpha)$-coalescent is proportional to
$N^{\alpha-1}$ generations ($1 < \alpha < 2$). The $N^{\alpha-1}$ unit
of time might require strong assumptions regarding the population size
or the mutation rate in order to recover the mutations in a given
sample used to estimate $\alpha$ when $\alpha$ is estimated to be near
1.





\bibliographystyle{alpha}
\begin{thebibliography}{99}





\bibitem[BLS15]{BLS15}
  \newblock \url{https://dx.doi.org/10.1214/18-ejp175}.
  \newblock {2018}.
\newblock {23}, {0}.
 \newblock {Matthias Birkner and Huili Liu and Anja Sturm}.
 \newblock {Coalescent results for diploid exchangeable population models},
   {Electronic Journal of Probability}.





\bibitem[D2024]{D2024}
 Diamantidis,  Dimitrios and Fan,  Wai-Tong (Louis) and Birkner,
 Matthias and Wakeley,  John.  {Bursts of coalescence within
 population pedigrees whenever big families occur}.  Genetics Volume 227,  February
  2024. \\
 \url{https://dx.doi.org/10.1093/genetics/iyae030}.



\bibitem[F2025]{F2025}
Frederic Alberti and Matthias Birkner and Wai-Tong Louis Fan and John
Wakeley. {A conditional coalescent for diploid exchangeable population
models given the pedigree}. {arXiv}, 2025
\\
\url{https://dx.doi.org/10.48550/arXiv.2505.15481}


\bibitem[CDEH25]{chetwyn-diggle_beta}
JA~Chetwyn-Diggle, Bjarki Eldon, and Matthias Hammer.
\newblock Beta-coalescents when sample size is large.
\newblock In preparation, 2025+.

\bibitem[DK99]{DK99}
P~Donnelly and T~G Kurtz.
\newblock Particle representations for measure-valued population models.
\newblock {\em Ann Probab}, 27:166--205, 1999.
\\
\url{https://dx.doi.org/10.1214/aop/1022677258}

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.


\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.


\bibitem[Pit99]{P99}
J~Pitman.
\newblock Coalescents with multiple collisions.
\newblock {\em Ann Probab}, 27:1870--1902, 1999.
\\
\url{https://dx.doi.org/10.1214/aop/1022874819}


\bibitem[Sag99]{S99}
S~Sagitov.
\newblock The general coalescent with asynchronous mergers of ancestral lines.
\newblock {\em J Appl Probab}, 36:1116--1125, 1999.
\\
\url{https://doi.org/10.1239/jap/1032374759}


\bibitem[Sch03]{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.
\\
\url{https://doi.org/10.1016/S0304-4149(03)00028-0}


\bibitem[S00]{S00}
J~Schweinsberg.
\newblock Coalescents with simultaneous multiple collisions.
\newblock {\em Electronic Journal of Probability}, 5:1--50, 2000.
\\
\url{https://dx.doi.org/10.1214/EJP.v5-68}



\bibitem[Tan11]{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.


\end{thebibliography}




@
\end{document}
