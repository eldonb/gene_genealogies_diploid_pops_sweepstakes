%%% copyright (c) 2024 bjarki eldon
%%% see annealederiNdiploid.cpp
%%% compare with quencheddiploidsfs.cpp
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
\usepackage{rotate}
\usepackage{orcidlink}
\setstretch{1}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\IN}{\ensuremath{\mathds N}}
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
%%
\newtheorem{thm}{Theorem}[section]
\newtheorem{propn}[thm]{Proposition}%
\newtheorem{lemma}[thm]{Lemma}%
\newtheorem{eg}[thm]{Example}
\newtheorem{defn}[thm]{Definition}
\newtheorem{remark}[thm]{Remark}
\newtheorem{notn}[thm]{Notation}
\newtheorem{corollary}[thm]{Corollary}
\newtheorem{conjecture}[thm]{Conjecture}
\newtheorem{assumption}[thm]{Assumption}
%%%%%%%%%

\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  \textsf{\textbf{\@@title}}
\end{flushleft}
\begin{center}
  \textsc{\@@author}
\end{center}
\egroup
}
\makeatother

\title{Gene genealogies in  diploid populations evolving according to sweepstakes reproduction \\ --- approximating $\EE{R_{i}^{N}(n)}$}
\author{Bjarki Eldon\footnote{\href{mailto:beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17
to Wolfgang Stephan; acknowledge funding by the Icelandic Centre of
Research (Rann\'is) through an Icelandic Research Fund
(Ranns\'oknasj\'o{\dh}ur) Grant of Excellence no.\ 185151-051 to Einar
\'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\ Etheridge, Wolfgang
Stephan, and BE. BE also acknowledges Start-up module grants through
SPP 1819 with Jere Koskela and Maite Wilke-Berenguer, and with Iulia
Dahmer. \\ \today} \orcidlink{0000-0001-9354-2391} }

\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


%%\rule{\textwidth}{.8pt}


\begin{abstract}
  Let $\mathds N_{0} \equiv \set{0,1,2, \ldots}$ and $\set{\xi^{n,N}}
\equiv \set{ \xi^{n,N}(g) : g \in \IN_{0} }$ be the ancestral process
tracking the random ancestral relations of a sample of gene copies;
$\tau^{N}(n) \equiv \inf \set{g \in \IN_{0} : \#\xi^{n,N}(g) = 1 } $
and $L_{i}^{N}(n) \equiv \sum_{j=0}^{\tau^{N}(n)} \# \set{ \xi \in
\xi^{n,N}(j) : \#\xi = i} $ for $i = 1,2, \ldots, 2n-1$, $L^{N}(n)
\equiv \sum_{j=0}^{\tau^{N}(n)} \# \xi^{n,N}(j) $, where $\xi^{n,N}(0)
\equiv \set{ \set{1,2}, \ldots, \set{2n-1,2n} } $. Then $R_{i}^{N}(n)
\equiv L_{i}^{N}(n)/L^{N}(n) $ and $L^{N}(n) =\sum_{j} L_{j}^{N}(n)$.
With this C++ code one estimates the functionals $\EE{R_{i}^{N}(n)}$
of gene genealogies of samples from a finite diploid panmictic
population of constant size absent selfing and evolving according to
specific models of sweepstakes reproduction
\end{abstract}



\tableofcontents




@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

This document and any source code it contains is distributed under the
terms of the GNU General Public Licence (version $\ge 3$).  You should
have received a copy of the licence along with this file (see file
COPYING).


The source codes described in this document are free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

This document and the code it contains is distributed in the hope that
it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.




@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \cite{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \cite{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library (GSL).
Compiles on {\tt Linux debian}  6.12.9-amd64  with {\tt ctangle}  4.11 and   {\tt g++} 14.2  and  {\tt GSL} 2.8

Using a Makefile can be helpful, naming this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        g++ -std=c++26 -O3 -march=native -m64 -x c++ iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:


{\tt valgrind -v ---leak-check=full ---leak-resolution=high ---num-callers=40 ---vgdb=full <program call>}


Use {\tt cppcheck} to check the code:

{\tt cppchek  ---enable=all ----language=c++ <prefix>.c}

To generate estimates on a computer with several CPUs it may be
convenient to put  in a text file ({\tt simfile}):


{\tt ./a.out \$(shuf -i 484433-83230401 -n1) > resout<i>}


for  $i = 1,\ldots, y$ and use  {\tt
parallel}\cite{tange22gnu_paral}

{\tt parallel ---gnu -jy :::: ./simfile}


@* {\bf intro}. 
\label{sec:intro}


Suppose a population is and evolves as in Definition~\ref{dschwpop}. 
\begin{defn}[Evolution of a diploid population]
\label{dschwpop}%
Consider a diploid (each diploid individual carries two gene copies or
chromosomes; each diploid individual can also be seen as a pair of
chromosomes) panmictic population of constant size $2N$ diploid
individuals.  In any given generation we arbitrarily label each
individual with a unique label, and form all possible $(2N-1)N$
(unordered) pairs of labels. Then we sample $N$ pairs of labels
independently and uniformly at random without replacement.  The $N$
parent pairs thus formed independently produce random numbers of
diploid potential offspring according to some given law.  Each
offspring receives two chromosomes, one chromosome from each of its
two parents, with each inherited chromosome sampled independently and
uniformly at random from among the two parent chromosomes.  If the
total number of potential offspring is at least $2N$, we sample $2N$
of them uniformly at random without replacement to survive to maturity
and replace the current individuals; otherwise we assume the
population is unchanged over the generation (all the potential
offspring perish before reaching maturity).
\end{defn}



\begin{remark}[Illustrating Definition~\ref{dschwpop}]
\label{rm:illdschwpop}
The mechanism described in Definition~\ref{dschwpop} is illustrated
below, where $\set{a,b}$ denotes a diploid individual carrying gene
copies $a,b$ and $\set{\set{a,b},\set{c,d}}$ denotes a pair of
parents. Here we have arbitrarily labelled the gene copies just for
the sake of illustrating the evolution over one generation.
\begin{displaymath}
\begin{split}
\text{stage} &\quad   \text{individuals involved} \\
1 & \quad  \set{\{a,b\} , \{c,d\}}, \ldots,  \set{\set{w,x}, \set{y,z}} \quad \text{:  $N$ parent pairs} \\
2 & \quad \underset{X_{1}}{\underbrace{\set{a,d},\set{b,d}, \ldots, \set{a,c}}}, \ldots, \underset{X_{N}}{ \underbrace{ \set{w,y}, \ldots, \set{x,z}} } \quad \text{ : $X_{1} + \cdots +  X_{N}$ potential offspring}  \\
3 & \quad    \set{b,d},\ldots, \set{x,y} \quad \text{ : $2N$ surviving offspring (whenever $X_{1}+ \cdots +X_{N} \ge 2N$)}
\end{split}
\end{displaymath}
In stage~1 above, the current $2N$ diploid individuals randomly form
$N$ pairs; in stage~2 the $N$ pairs formed in stage~1 independently
produce random numbers $X_{1}, \ldots, X_{N}$ of potential offspring,
where each offspring receives one gene copy (or chromosome) from each
of its two parents; in the third stage $2N$ of the
$X_{1} + \cdots + X_{N}$ potential offspring (conditional on there
being at least $2N$ of them) are sampled uniformly and without
replacement to survive to maturity and replace the parents.
\end{remark}
It follows from Definition~\ref{dschwpop}, and emphasized in
Remark~\ref{rm:illdschwpop}, that two gene copies in the same diploid
individual must have come from separate parents . It is thus necessary
to track the pairing of marked gene copies (ancestral to the sampled
gene copies) in diploid individuals.


Let $\set{\xi^{n,N}(r) : r\in \mathds N\cup\set{0} }$ denote the
{\it ancestral process}, a Markov sequence tracking the ancestral
relations of sampled gene copies in a finite population. Let
$\set{\xi^{n}(t) : t \ge 0}$ denote a coalescent, a Markov chain on
the partitions of $[2n]$.
Let $ \# A$ be the cardinality of a given set $A$, write 
$\tau^{N}(n) := \inf \set{j \in \mathds N : \# \xi^{n,N}(j) = 1}$, and 
$\tau(n) := \inf \set{t \ge 0 : \# \xi^{n}(t) = 1 }$. Consider the functionals  
\begin{displaymath}
\begin{split}
L_{i}^{N}(n) := \sum_{j=1}^{\tau^{N}(n)}  \# \set{\xi \in \xi^{n,N}(j) : \# \xi   = i }, &  \quad   L^{N}(n) :=  \sum_{j=1}^{\tau^{N}(n)} \# \xi^{n,N}(j) \\
L_{i}(n) := \int_{0}^{\tau(n)} \# \set{\xi \in \xi^{n}(t) : \# \xi = i } dt, &  \quad L(n) :=  \int_{0}^{\tau(n)} \# \xi^{n}(t) dt
\end{split}
\end{displaymath}
\cite{Birkner2018}.  The functionals $L_{i}^{N}(n)$ and $L_{i}(n)$ are the
random length of branches supporting $i\in [n-1]$ leaves,
$L^{N}(n) = L_{1}^{N}(n) + \cdots + L_{n-1}^{N}(n)$, and
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$.




Recall $n$ is the number of diploid individuals sampled, so we have
$2n$ gene copies.  Let $L_{i}^{N}(n)$ denote the random length of
branches supporting $i \in \{1,2,\ldots, 2n-1\}$ leaves; write
$R_{i}^{N}(n) \equiv L_{i}^{N}(n)/\sum_{j=1}^{2n-1}(n)$.  We estimate
$\EE{R_{i}^{N}(n)}$ when the sample comes from a finite haploid
panmictic population evolving according to random recruitment
success. By ``recruitment success'' we mean an increased chance of
producing surviving offspring that number on the order of the
population size.




Write $n = \{1,2, \ldots, n\}$ for any
$n \in \mathds{N} := \{1,2,\ldots \}$.  Let $X_{1}, \ldots, X_{N}$
denote the random number of potential offspring produced by the
current individuals. They are independent but may not always be
identically distributed.  Let $X$ be independent of
$X_{1}, \ldots, X_{N}$ and suppose
\begin{equation}
\label{eq:pmf}
\prb{X=k} =  C\left( \frac{1}{k^{a}} -  \frac{1}{(1+k)^{a}}  \right), \quad k \in \set{1,2,\ldots, \zeta(N)}
\end{equation}
with $a > 0$ and where $C > 0$ is such that  $\prb{X \in [\zeta(N)]} = 1$.  Write
$X\vartriangleright L(a,\zeta(N))$ when $X$ is distributed according
to \eqref{eq:pmf} for some given $a$ and $\zeta(N)$.  In any given
generation the current individuals independently produce potential
offspring according to \eqref{eq:pmf}; from the pool of potential
offspring $N$ of them are sampled uniformly at random and without
replacement to survive and reach maturity and replace the current
individuals.  The model in \eqref{eq:pmf} is an extension of the one in \cite[Equation~11]{schweinsberg03}.

Let $0 < \alpha < 2$ and $\kappa \ge 2$ be fixed.  We consider two
models.  Let $E$ be the event 
\begin{displaymath}
E \equiv \set{X_{i} \vartriangleright \mathds L(\alpha,\zeta(N)), i = 1,2,\ldots, N }
\end{displaymath}
and $E^{\sf c}$ where $\kappa$ replaces  $\alpha$ in $E$.  
The  $X_{1}, \ldots, X_{N}$ are iid and 
\begin{equation}
\label{eq:model0}
X_{1} \vartriangleright \mathds L\left( \one{E } \alpha +   \one{E^{\sf c} } \kappa, \zeta(N)\right) 
\end{equation}


We also consider another model where the $X_{1}, \ldots, X_{N}$ are
independent, $X_{i}$ is as in \eqref{eq:model0} for $i$ picked
uniformly at random and
$X_{j\neq i} \vartriangleright L(\kappa,\zeta(N))$.  In this second
model the $X_{1}, \ldots, X_{N}$ are independent but may not always be
identically distributed. For suitable choices of $\varepsilon_{N}$ one
can identify the limiting diffusions/coalescents. Clearly, choosing
$\varepsilon_{N} \ge 1$ results in iid $X_{1},\ldots X_{N}$ with
$X_{1} \vartriangleright L(\alpha, \zeta(N))$; taking
$\varepsilon_{N} \le 0$ gives an alternative to the Wright-Fisher
model.


The algorithm is summarized in \S~\ref{sec:code}; the code follows in 
\S~\ref{sec:includes}--\S~\ref{sec:main}; we conclude in
\S~\ref{sec:concl}

@* {\bf code}.
\label{sec:code}

We sample $n$ diploid individuals and so $2n$ gene copies; the
algorithm summarized :

\begin{enumerate}
\item $\svigi{r_{1}^{N}(n), \ldots, r_{2n-1}^{N}(n)} \leftarrow \svigi{0,\ldots, 0}$
\item For each of $M$ experiments \S~\ref{sec:one}
\begin{enumerate}
\item $\left(\ell_{1}^{N}(n), \ldots, \ell_{n-1}^{N}(n) \right) \leftarrow (0,\ldots, 0)$
\item initialise block sizes  $\svigi{ ({b_{1},b_{2}}), \ldots, (b_{2n-1}, b_{2n}) } \leftarrow  \svigi{ \svigi{1,1}, \ldots, \svigi{ 1,1} } $
\item $m\leftarrow 2n$ the current number of ancestral gene copies
\item {\bf while} $m > 1$ :
\begin{enumerate}
\item $\ell_{b_{j}}^{N}(n) \leftarrow 1 +  \ell_{b_{j}}^{N}(n)  $ for $b_{j}=b_{1},\ldots, b_{m}$ \S~\ref{sec:updatebi}
\item sample numbers of potential offspring $X_{1},\ldots, X_{N}$ \S~\ref{sec:sampleSN} 
\item given $X_{1},\ldots, X_{N}$  assign marked diploid individuals to families  \S~\ref{sec:rmvhyper} 
\item merge blocks and append new blocks to sample
\S~\ref{sec:mergeblocks}
\end{enumerate}
\item given branch lengths $\ell_{i}^{N}(n)$ update estimate of  $\EE{R_{i}^{N}(n)}$ \S~\ref{sec:updateri}
\begin{displaymath}
r_{j}^{N}(n) \leftarrow  \frac{\ell_{j}^{N}(n) }{ \sum_{i}\ell_{i}^{N}(n) } +  r_{j}^{N}(n), \quad j = 1, \ldots, 2n-1
\end{displaymath}
\end{enumerate}
\item return an  estimate $(1/M) r_{i}^{N}(n)$ of   $\EE{R_{i}^{N}(n)}$ \S~\ref{sec:estimate} 
\end{enumerate}


@*1 {\bf includes}.
\label{sec:includes}

the included libraries; we use the GSL library so this needs to be
installed

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
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <unordered_set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>


@*1 {\bf constants}.
\label{sec:constants}


the parameter values

@<constants@>=@#
const double CONST_ALPHA = 1.01 ;
const double CONST_A = 2.0 ;
/* *** \newline  $N$ is  number of parent pairs producing potential offspring;
   $2N$  is number of diploid individuals
   and so $4N$  gene copies in the population
   *** */
const std::size_t CONST_N = 1e3 ;
const double dN = static_cast<double>(CONST_N); 
const std::size_t CONST_CUTOFF = 2*CONST_N*static_cast< std::size_t>( ceil( log(dN)) ) ;
const std::vector<double> P =  {0.25, 0.5, 0.75, 1.};
const double CONST_VAREPSILON = 0.1 ;
/*  \newline  initialise number of diploid individuals sampled ;
   $n$ is number of diploid individuals sampled so
   $2n$  is number of gene copies (blocks) *** */
const std::size_t CONST_SAMPLE_SIZE = 1e3 ;
const int CONST_NUMBER_EXPERIMENTS = 2000 ;


@*1 {\bf the random number generators}.
\label{sec:rngs}

the GSL random number generator 


@<rngs@>=@#
gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the probability mass function for  potential offspring}.
\label{sec:pmf}

the probability mass function \eqref{eq:pmf}

@<pmf@>=@#
static double kernel( const double &k,   const double &a  )
{
  /* \newline compute $k^{-a} - (1+k)^{-a}$ */
  return ( pow( 1./k, a) - pow( 1./(1. + k), a) ) ;
}


@*1 {\bf generate cmfs}.
\label{sec:cmfs}

generate cumulative mass functions 

@<cmfs@>=@#
static void generatecmf( std::vector<double>& vcmfalpha, std::vector<double>& vcmfA )
{
  /* \newline  the cmf arrays |vcmfalpha| and |vcmfA|  are initialised to  $(0, ..., 0)$  *** */
  double salpha {} ;
  double sA {} ;
  for( std::size_t i = 2 ; i < CONST_CUTOFF; ++i){
  /* \newline |kernel| \S~\ref{sec:pmf} */
    vcmfalpha[i] = vcmfalpha[i-1] + kernel( static_cast<double>(i), CONST_ALPHA);
    salpha +=  kernel( static_cast<double>(i), CONST_ALPHA);
    vcmfA[i] = vcmfA[i-1] + kernel( static_cast<double>(i), CONST_A);
    sA +=  kernel( static_cast<double>(i), CONST_A);}
  assert(salpha > 0.);
  assert( sA > 0.);
  std::transform( vcmfalpha.begin(), vcmfalpha.end(), vcmfalpha.begin(), [&salpha]( const auto&x ){ return x/salpha; });
  std::transform( vcmfA.begin(), vcmfA.end(), vcmfA.begin(), [&sA]( const auto&x ){ return x/sA; });

  vcmfalpha[ CONST_CUTOFF] = 1. ;
  vcmfA[CONST_CUTOFF] = 1. ;
}



@*1 {\bf sample one  random number of potential offspring}.
\label{sec:randomX}

sample one random number of potential offspring using \eqref{eq:pmf};
return $\min\{j : u \le F(j) \}$ where $F$ is a cumulative mass function and $u$ a random uniform 

@<sample random $X$@>=@#
static std::size_t randomX( const std::vector<double>& f)
{
  const double u = gsl_rng_uniform(rngtype);
  std::size_t j {2} ;
  while( u > f[j]){ ++j ; }

  return j ;
}


@*1 {\bf sample a pool of potential offspring}.
\label{sec:sampleSN}

sample a pool of potential offspring using \S~\ref{sec:randomX}

@<sample pool potentials@>=@#
static std::size_t sampleSN(  std::vector<std::size_t>& vX, const std::vector<double>& f )
{
  /* \newline  sample a pool  of potential offspring */
  std::size_t sn {} ; 
  for( std::size_t i = 0 ; i < CONST_N ; ++i){
  /* \newline |randomX| \S~\ref{sec:randomX} */
     vX[i] = randomX(f);
    
    sn += vX[i] ;}

  return sn ;
}



@*1 {\bf sample number per family}.
\label{sec:rmvhyper}

Given a pool of potential offspring sample a random number of
surviving marked diploid offspring assigned to a family; a `marked'
diploid offspring is one carrying a gene copy (chromosome) ancestral
to a sampled chromosome 

@<random multivariate hypergeometric@>=@#
static void rmvhyper( std::vector<std::size_t>& merger_sizes, std::size_t&  k, const std::vector<size_t>& v_juvs, const std::size_t& SN )
{
  /* \newline  |k| is the current number of marked individuals  */
  merger_sizes.clear();
  
  size_t n_others =  SN - v_juvs[0] ;
  /* \newline  sample the number of marked diploid potential offspring  assigned to the first family */
  size_t x = gsl_ran_hypergeometric( rngtype, v_juvs[0], n_others, k);
  if( x > 0){
    /* \newline  only record  merger sizes */
    merger_sizes.push_back(x ); }
    /* \newline  update the remaining number of blocks  */
  k -= x ;
  /* \newline  update new number of lines */
  
  size_t i =0 ;
    /* \newline  we can stop as soon as all lines 
       have been assigned to a family */
  while( (k > 0) && (i < CONST_N - 1) ){
  /* \newline  set the index to the one being sampled from  */
    ++i ;
    /* \newline  update |n_others|  */
    n_others -= v_juvs[i] ;
    x =  gsl_ran_hypergeometric( rngtype, v_juvs[i], n_others, k );
    if( x > 0){
      merger_sizes.push_back( x) ; }
      /* \newline  update the remaining number of blocks  */
      k -= x ;
  }
  /* \newline  check if at least two lines assigned to last individual  */
  if( k > 0){
    merger_sizes.push_back( k); }
}



@*1 {\bf sample parent chromosome}.
\label{sec:samplepchrom}

assigna a block to one of four  parent chromosomes

@<sample parent chromosome@>=@#
static std::size_t samplepchrom()
{
  const double u = gsl_rng_uniform(rngtype);

  return (u <0.25 ? 0 : (u < .5 ? 1 : (u < 0.75 ? 2 : 3) ) ) ;  
}


@*1 {\bf merge blocks}.
\label{sec:mergeblocks}

merge blocks; blocks in the same diploid individual cannot merge, they
disperse and have to be assigned to separate parents

@<merge blocks@>=@#
static void mergeblocks( const std::size_t& j, std::vector< std::pair< std::size_t, std::size_t>>& sample, const std::size_t &x)
{
  /* \newline 
     |sample| is a vector of pairs of block sizes; need to keep track of pairing of
     blocks in diploid individuals; blocks in same diploid individual disperse and
     so cannot merge
  */

  /* \newline  |x| is number of marked diploid  offspring of a  given family */

  /* \newline |newblocks| will store the new block sizes of merged and/or continuing blocks */
  std::vector<std::size_t> newblocks (4, 0);
  assert( x > 0);

 /* \newline |sample| is a vector of marked diploid  individuals  */

  assert( j+x <= sample.size());
  
  for( std::size_t i = j ; i < j+x; ++i){
  /* \newline if marked diploid  individual  |i| is double-marked then the blocks disperse; otherwise use |samplepchrom| \S~\ref{sec:samplepchrom} to sample parent chromosome */
    newblocks[ sample[i].second > 0 ? ( gsl_rng_uniform(rngtype) < 0.5 ? 0 : 1) : samplepchrom() ] += sample[i].first ;
    newblocks[ sample[i].first > 0 ? ( gsl_rng_uniform(rngtype) < 0.5 ? 2 : 3) : samplepchrom() ] += sample[i].second ;}



 /* \newline append  the new blocks in pairs  to the sample */

  if( (newblocks[0] + newblocks[1]) > 0 ){
    sample.push_back( std::make_pair( newblocks[0], newblocks[1])) ; }
  
  if( (newblocks[2] + newblocks[3]) > 0 ){
    sample.push_back( std::make_pair( newblocks[2], newblocks[3])) ; }
}


@*1 {\bf update sample}.
\label{sec:updatesample}


update sample by merging blocks and appending  the new blocks to the sample 

@<update sample@>=@#
static std::size_t  updatesample( const std::size_t &offset, const std::vector< std::size_t>& vnu, std::vector<std::pair< std::size_t, std::size_t>>& sample )
{
  /* \newline   |vnu|  is number of marked diploid offspring  per parent pair */

  const std::size_t newoffset = sample.size() ;

  /* \newline  |o| is the index of the first element to be worked with */
  std::size_t o = offset ; 
  for( const auto &m : vnu){
    assert( m > 0);
    /* \newline \S~\ref{sec:mergeblocks} */
    mergeblocks( o, sample, m);
    o += m ;}

/* \newline |newoffset| is the index of the first element of the new sample */
  return newoffset ;
}


@*1 {\bf update site-frequency spectrum $\ell_{i}^{N}(n)$}. 
\label{sec:updatebi}

update the site-frequency spectrum (branch lengths $\ell_{i}^{N}(n)$);
the sample is a vector of pairs of block sizes 

@<update branch lengths $\ell_{i}^{N}(n)$@>=@#
static void updatebi( std::vector<std::size_t>& b, const std::vector<std::pair< std::size_t, std::size_t>>& s )
{
  for( const auto &p: s){
  assert( p.first < 2*CONST_SAMPLE_SIZE );
    b[0] += (p.first > 0 ? 1 : 0);
    b[p.first] += (p.first > 0 ? 1 : 0);
    assert(p.second < 2*CONST_SAMPLE_SIZE );
    b[0] += (p.second > 0 ? 1 : 0);
    b[p.second] += (p.second > 0 ? 1 : 0);}
}



@*1 {\bf update estimate $r_{i}^{N}(n)$ of $\EE{R_{i}^{N}(n)}$}.
\label{sec:updateri}

given a realisation of branch lengths $\ell_{i}^{N}(n)$ update
estimate $r_{i}^{N}(n)$ of $\EE{R_{i}^{N}(n)}$ 

@<update $r_{i}^{N}(n)$@>=@#
static void updateri( std::vector<double>& ri, const  std::vector<std::size_t>& vb )
{
  assert( vb[0] > 0);
  for( std::size_t i = 1; i < (2*CONST_SAMPLE_SIZE); ++i){
    ri[i] += (static_cast<double>( vb[i])/static_cast<double>( vb[0]) ) ; }
}




@*1 {\bf one realisation of branch lengths}.
\label{sec:one}

generate one realisation of branch lengths $\ell_{i}^{N}(n)$


@<one branch length realisation@>=@#
static void one( const std::vector<double>& vfalpha, const std::vector<double>& vfA, std::vector< std::size_t>& vbi, std::vector<double>& vri )
{

  /* \newline  initialise sample $((1,1),\ldots, (1,1))$ */
  std::vector< std::pair< std::size_t, std::size_t>> sample (CONST_SAMPLE_SIZE, std::make_pair(1,1) );

  std::vector< std::size_t> vx (CONST_N, 0);
  std::vector< std::size_t> vm {} ;
  vm.reserve( 2*CONST_SAMPLE_SIZE);
  std::size_t npo {} ;
  std::size_t ss {} ;
  std::size_t offs = 0 ;
  std::fill( std::begin(vbi), std::end(vbi), 0);
  
   while( (sample.back().first < 2*CONST_SAMPLE_SIZE) && (sample.back().second < 2*CONST_SAMPLE_SIZE))
    {
    /*  \newline  |updatebi|  \S~\ref{sec:updatebi} */
      updatebi( vbi, sample);

      /* \newline  sample potential offspring  |sampleSN|  \S~\ref{sec:sampleSN}  */
      npo = sampleSN( vx, (gsl_rng_uniform(rngtype) < CONST_VAREPSILON ? vfalpha : vfA) );

      /*  \newline  sample number of marked potential offspring per family */
      /* \newline |ss| is number of marked diploid individuals */
      ss = sample.size() - offs ;
      /* \newline |rmvhyper|  \S~\ref{sec:rmvhyper} */
      rmvhyper( vm, ss, vx, npo) ;
      assert( static_cast<std::size_t>(std::accumulate( vm.begin(), vm.end(),0)) == (sample.size() - offs) ) ;
      /* \newline merge blocks and append new blocks to sample  |updatesample|  \S~\ref{sec:updatesample} */
      offs = updatesample(offs, vm, sample);
    }

/* \newline given realisation of branch lengths update estimate $r_{i}^{N}(n)$ \S~\ref{sec:updateri} */
   updateri( vri, vbi);
}



@*1 {\bf estimate $\EE{R_{i}^{N}(n)}$}.
\label{sec:estimate}

get estimates $r_{i}^{N}(n)$ of  $\EE{R_{i}^{N}(n)}$


@<estimate $\EE{R_{i}^{N}(n)}$@>=@#
static void estimate()
{
  std::vector<double> cmfalpha (CONST_CUTOFF + 1, 0.);
  std::vector<double> cmfA (CONST_CUTOFF + 1, 0.);
  /* \newline |generatecmf| \S~\ref{sec:cmfs} */
  generatecmf(cmfalpha, cmfA);

  std::vector< std::size_t> branchlengths (2*CONST_SAMPLE_SIZE, 0);
  std::vector< double> eri (2*CONST_SAMPLE_SIZE, 0.);

  int r = CONST_NUMBER_EXPERIMENTS + 1 ;
  while( --r > 0){
  /* \newline |one|  \S~\ref{sec:one} */
    one( cmfalpha, cmfA, branchlengths, eri ) ;}

  for( const auto &e : eri){ std::cout << e << '\n';}
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
/* \newline \S~\ref{sec:pmf} */
@<pmf@>@#
/* \newline \S~\ref{sec:cmfs} */
@<cmfs@>@#
/* \newline \S~\ref{sec:randomX} */
@<sample random $X$@>@#
/* \newline \S~\ref{sec:sampleSN} */
@<sample pool potentials@>@#
/* \newline \S~\ref{sec:rmvhyper} */
@<random multivariate hypergeometric@>@#
/* \newline \S~\ref{sec:samplepchrom} */
@<sample parent chromosome@>@#
/* \newline \S~\ref{sec:mergeblocks} */
@<merge blocks@>@#
/* \newline \S~\ref{sec:updatesample} */
@<update sample@>@#
/* \newline \S~\ref{sec:updatebi} */
@<update branch lengths $\ell_{i}^{N}(n)$@>@#
/* \newline \S~\ref{sec:updateri} */
@<update $r_{i}^{N}(n)$@>@#
/* \newline \S~\ref{sec:one} */
@<one branch length realisation@>@#
/* \newline \S~\ref{sec:estimate} */
@<estimate $\EE{R_{i}^{N}(n)}$@>@#



int main(int argc, const char * argv[])
{

  /* \newline \S~\ref{sec:rngs} */
  setup_rng( static_cast< std::size_t>(atoi(argv[1])) ) ;

/* \newline \S~\ref{sec:estimate} */
  estimate();

  gsl_rng_free(rngtype);
  return GSL_SUCCESS; 
}



@* {\bf conclusions and bibliography}.
\label{sec:concl}


We estimate mean relative branch lengths $\EE{R_{i}^{N}(n)}$ for a
sample from a finite diploid panmictic population of constant size
(see Definition~\ref{dschwpop} and Remark~\ref{rm:illdschwpop})  and
evolving according to random recruitment success \eqref{eq:pmf}.   Since the
population is diploid each individual can carry two ancestral blocks;
we exclude selfing (clonal reproduction) so blocks in the same diploid
individual cannot merge; they disperse to the two parents.


One would want to compare $\EE{R_{i}^{N}(n)}$ for some pre-limiting
model to $\EE{R_{i}(n)}$ predicted by the coalescent for which the
pre-limiting model is in the domain-of-attraction to, in an effort to
investigate how well the trees predicted by the two processes agree.





\begin{thebibliography}{BLS18}

\bibitem[BLS18]{Birkner2018}
Matthias Birkner, Huili Liu, and Anja Sturm.
\newblock Coalescent results for diploid exchangeable population models.
\newblock {\em Electronic Journal of Probability}, 23(none), January 2018.

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.

\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.

\bibitem[Sch03]{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.

\bibitem[Tan22]{tange22gnu_paral}
O~Tange.
\newblock {\em GNU {Parallel} 20221122}, 2022.
\newblock Zenodo.

\end{thebibliography}




@
\end{document}
