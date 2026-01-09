%%% copyright (c) 2024 bjarki eldon
%%% see annealederiNdiploid.cpp
%%% compare with quencheddiploidsfs.cpp
\pdfoutput=1
\documentclass[a4paper,10pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
%%\usepackage[lf]{Baskervaldx}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
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
\usepackage{orcidlink}
\setstretch{1}
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\small \textsf{\bf  #1 }}}
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

\title{Gene genealogies for diploid populations evolving according to sweepstakes reproduction  \\ --- approximating $\EE{\widetilde R_{i}^{N}(n)}$}
\author{Bjarki Eldon\footnote{ \href{beldon11@@gmail.com}{beldon11@@gmail.com}} \footnote{Supported by 
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
With this C++ code one estimates mean conditional relative branch
lengths $ \rho_{i}^{N}(n) \equiv \EE{\widetilde R_{i}^{N}(n)}$ where
$R_{i}^{N}(n) = L_{i}^{N}(n)/\sum_{j=1}^{2n-1}L_{j}^{N}(n)$ and
$L_{i}^{N}(n)$ is the random length of branches supporting
$i\in \{1,2,\ldots, 2n-1\}$ leaves where $n$ is the number of diploid
individuals sampled (hence $2n$ gene copies); $\mathcal{A}^{(N,n)}$ is
the sigma field keeping track of the ancestral relations of whole
population.  The idea is that the gene copies of a sample share an
ancestry, a gene genealogy, even if the genealogy is unknown. The
quantity $ \EE{\widetilde R_{i}^{N}(n)} $ is mean relative branch
length conditional on the population ancestry, the ancestral relations
of all gene copies at all times.  The sample comes from a finite
diploid panmictic population of constant size evolving according to
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
.c} file), one needs the GNU Scientific Library.

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
\tab\quad\quad\quad\quad        g++ -std=c++26 -O3 -march=native -m64 -xc++  iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:


{\tt valgrind -v ---leak-check=full ---leak-resolution=high ---num-callers=40 ---vgdb=full <program call>}


Use {\tt cppcheck} to check the code:

{\tt cppchek  ---enable=all ----language=c++ <prefix>.c}

To generate estimates on a computer with several CPUs it may be
convenient to put  in a text file ({\tt simfile}):


{\tt ./a.out \$(shuf -i <smallest random seed>-<largest random seed>  -n1) > resout<i>}

for  $i = 1,\ldots, y$ for suitables choices of random seeds so that no two replications  will have the same seed     and use  {\tt
parallel}\cite{tange22gnu_paral}

{\tt parallel ---gnu -jy :::: ./simfile}


@* {\bf intro}.
\label{sec:intro}



Suppose a population is and evolves as in Definition~\ref{dschwpop}. 
\begin{defn}[Evolution of a diploid population]
\label{dschwpop}
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



Write $[n] \equiv \set{1,2,\ldots, n}$ for $n\in \mathds N \equiv \set{1,2, \ldots}$.  
Let $\set{\xi^{n,N}(r) : r\in \mathds N\cup\set{0} }$ denote the
{\it ancestral process}, a Markov sequence tracking the ancestral
relations of sampled gene copies in a finite population. Let
$\set{\xi^{n}(t) : t \ge 0}$ denote a coalescent, a Markov chain on
the partitions of $[2n]$.
Let $ \# A$ be the cardinality of a given set $A$, write 
$\tau^{N}(n) := \inf \set{j \in \mathds N : \# \xi^{n,N}(j) = 1}$, and 
$\tau(n) := \inf \set{t \ge 0 : \# \xi^{n}(t) = 1 }$. Consider the functionals, for $i \in [2n-1]$,   
\begin{displaymath}
\begin{split}
L_{i}^{N}(n) := \sum_{j=1}^{\tau^{N}(n)}  \# \set{\xi \in \xi^{n,N}(j) : \# \xi   = i }, &  \quad   L^{N}(n) :=  \sum_{j=1}^{\tau^{N}(n)} \# \xi^{n,N}(j) \\
L_{i}(n) := \int_{0}^{\tau(n)} \# \set{\xi \in \xi^{n}(t) : \# \xi = i } dt, &  \quad L(n) :=  \int_{0}^{\tau(n)} \# \xi^{n}(t) dt
\end{split}
\end{displaymath}
\cite{Birkner2018}.  The functionals $L_{i}^{N}(n)$ and $L_{i}(n)$ are
the random length of branches supporting $i\in [2n-1]$ leaves,
$L^{N}(n) = L_{1}^{N}(n) + \cdots + L_{n-1}^{N}(n)$, and
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$.



Write $R_{i}^{N}(n) \equiv
L_{i}^{N}(n)/\sum_{j=1}^{2n-1}L_{j}^{N}(n)$.  We estimate
$\EE{\widetilde R_{i}^{N}(n)}$ when the sample comes from a finite
haploid panmictic population evolving according to sweepstakes
reproduction (heavy-tailed offspring number distribution).



Let $X_{1}, \ldots, X_{N}$ denote the random number of potential
offspring produced by the current individuals. They are independent
but may not always be identically distributed.  Let $X$ be independent
of $X_{1}, \ldots, X_{N}$ and suppose
\begin{equation}
\label{eq:pmf}
\prb{X=k} =  C\left( \frac{1}{k^{a}} -  \frac{1}{(1+k)^{a}}  \right), \quad k \in \set{2,3, \ldots, \zeta(N)}
\end{equation}
with $a > 0$ and where $\prb{X \in \{2,3, \ldots, \zeta(N) \} } = 1$.
Write $X\vartriangleright L(a,\zeta(N))$ when $X$ is distributed
according to \eqref{eq:pmf} for some given $a$ and $\zeta(N)$.  In any
given generation the current individuals independently produce
potential offspring according to \eqref{eq:pmf}; from the pool of
potential offspring $N$ of them are sampled uniformly at random and
without replacement to survive and reach maturity and replace the
current individuals. Note that \eqref{eq:pmf} is an extension of   \cite[Equation~11]{schweinsberg03}.


Let $0 < \alpha < 2$ and $\kappa \ge 2$ be fixed.  We consider two
models; one where $X_{1}, \ldots, X_{N}$ are iid and in generation $g$
with $(U_{g})_{g}$ a sequence of iid uniforms
\begin{equation}
\label{eq:model0}
X_{1} \vartriangleright L\left( \one{U_{g} \le \varepsilon_{N}} \alpha +   \one{U_{g} > \varepsilon_{N}} \kappa, \zeta(N)\right)
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


The quantity $\EE{\widetilde R_{i}^{N}(n)}$ is the mean relative
branch length conditioned on the population ancestry (the ancestral
relations of all the gene copies in the population at all times).  The
idea is that the gene copies in a sample are related through a single
tree, a gene genealogy, even if the tree is unknown.  We approximate
$\EE{\widetilde R_{i}^{N}(n)}$ by averaging over population
ancestries.

The algorithm is summarized in \S~\ref{sec:code}, the code follows in
\S~\ref{sec:includes}--\S~\ref{sec:main}, we conclude in  \S~\ref{sec:concl}. Comments within the code are in \aths{this font and color}


@* {\bf code}.
\label{sec:code}


The ancestry $\mathds{A} = (A_{i}(g))_{i,g}$ records the ancestral
relations of the gene copies (chromosomes); if the chromosome on level
$\ell$ at time $g$ produces $k$ surviving copies then
$A_{j_{1}}(g+1) = \ldots = A_{j_{k}}(g+1) = \ell$; and if
$A_{\ell}(g) = A_{k}(g)$ then the individuals living on level $\ell$
resp.\ $k$ at time $g$ share an immediate ancestor.  The ancestry
$\mathds{A}$ records the ancestral relations of $4N$ gene copies for
each generation, where $N$ is the number of parent pairs producing
potential offspring.



\begin{enumerate}
\item $\svigi{\widehat r_{i}^{N}(n), \ldots, \widehat r_{2n-1}^{N}(n)   }\leftarrow (0, \ldots ,0) $
\item For each of $M$ experiments \S~\ref{sec:estimate} :
\begin{enumerate}
\item initialise the ancestry  $A_{\ell}(0) = \ell$ for $\ell = 1,2, \ldots,  4N$ 
\item $\left( \ell_{1}^{N}(n), \ldots, \ell_{n-1}^{N}(n) \right) \leftarrow (0, \ldots, 0)$
\item {\bf while} a sample tree is incomplete \S~\ref{sec:getcompletesampletree} : 
\begin{enumerate}
\item add to the ancestry $\mathds{A}$  \S~\ref{sec:addgeneration} by sampling a pool of potential offspring using \S~\ref{sec:randomX} and assigning gene copies (chromosomes) to each offspring \S~\ref{sec:assignchromstooffspring} 
\item take a random sample \S~\ref{sec:randomsample} 
\item check if the sample is complete \S~\ref{sec:treeincomplete} 
\end{enumerate}
\item given a complete sample tree read the branch lengths $\ell_{i}^{N}(n)$ off the fixed tree  \S~\ref{sec:readcompletetree} 
\item update estimate $\widehat r_{i}^{N}(n)$   of     $\EE{ \widetilde R_{i}^{N}(n) }$  \S~\ref{sec:updaterrs}
\begin{displaymath}
\widehat r_{i}^{N}(n) \leftarrow  \frac{\ell_{i}^{N}(n) }{ \sum_{j} \ell_{j}^{N}(n) } + \widehat r_{i}^{N}(n)
\end{displaymath}
\end{enumerate}
\item return the  estimates   $(1/M)\widehat r_{i}^{N}(n)$ of   $ \EE{ \widetilde R_{i}^{N}(n) } $
\end{enumerate}


@*1 {\bf includes}.
\label{sec:includes}

the included libraries; we use the GSL library

@<includes@>=@#
#include <iostream>
#include <fstream>
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
#include <chrono>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <unordered_set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>


@*1 {\bf constants}.
\label{sec:constants}


the parameter values

@<constants@>=@#
/* \newline \aths{ $N$ diploid individuals so $2N$ gene copies} */
const std::size_t CONST_N = 1e2; 
const std::size_t CONST_POP_SIZE = 2*CONST_N;
const double dN = static_cast<double>( CONST_N);
/* \newline \aths{the upper bound $\zeta(N)$} */
const std::size_t CONST_CUTOFF = CONST_POP_SIZE * static_cast<std::size_t>( ceil(log(dN))) ;
/* \newline \aths{ $\alpha$} */
const double CONST_ALPHA = 1.01;
/* \newline \aths{ $\kappa$} */
const double CONST_A = 2.0 ;
/* \newline \aths{ $\varepsilon_{N}$ the probability of favorable conditions} */
const double CONST_VAREPSILON = 0.1 ;
/*  \newline \aths{ sample size $n$ of diploid individuals then $2n$ gene copies in sample} */
const std::size_t CONST_SAMPLE_SIZE = 1e1 ;
const int CONST_EXPERIMENTS = 1 ;



@*1 {\bf the random number generators}.
\label{sec:rngs}

the STL and GSL random number generators 


@<rngs@>=@#
/* \newline  \aths{ obtain a seed out of thin air for the random number engine} */
  std::random_device randomseed;
  /* \newline \aths{  Standard Mersenne twister  random number engine seeded with |randomseed()| } */
  std::mt19937_64 rng(randomseed());

gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the probability  mass function for number of offspring}.
\label{sec:pmf}

the probability mass function \eqref{eq:pmf} 

@<pmf@>=@#
static double  pmf( const std::size_t &k, const double& a )
{
  /* \newline \textcolor{violet}{\small \sf return $k^{-a} - (1+k)^{-a}$} */
  return ( pow(1./static_cast<double>(k), a) - pow( 1./static_cast<double>(k+1), a) ) ;
}


@*1 {\bf generate a cumulative mass function}.
\label{sec:generatecmf}

generate a cumulative mass function for sampling numbers of potential offspring


@<generate cmf@>=@#
static void generatecmf( std::vector<double>& cmfalpha,  std::vector<double> & cmfA)
{
  /* \newline \textcolor{violet}{\small \sf   the containers are initialised to $(0, ..., 0)$ }  */
  double salpha {} ;
  double sA {} ;
  for( std::size_t i = 2; i <= CONST_CUTOFF; ++i){
  /* \newline  \textcolor{violet}{\small \sf  |pmf|   \S~\ref{sec:pmf} } */
    cmfalpha[i] = cmfalpha[i-1] + pmf(i, CONST_ALPHA);
    cmfA[i] = cmfA[i-1] + pmf(i, CONST_A);
    salpha += pmf(i, CONST_ALPHA);
    sA += pmf(i, CONST_A);}

  assert(salpha > 0.);
  assert(sA > 0.);

  std::transform( cmfalpha.begin(), cmfalpha.end(), cmfalpha.begin(), [&salpha](const auto &x){return x/salpha;} );
  std::transform( cmfA.begin(), cmfA.end(), cmfA.begin(), [&sA](const auto &x){return x/sA;} );
  cmfalpha[CONST_CUTOFF] = 1. ;
  cmfA[CONST_CUTOFF] = 1. ;
}



@*1 {\bf sample a random number of potential offspring}.
\label{sec:randomX}

sample one random number of potential offspring
$\min\{ j : u \le F(j)\}$ where $u$ a random uniform and $F$ the
cumulative mass function


@<get one random number of potential offspring@>=@#
static std::size_t randomX( const  std::vector<double>& f )
{
  const double u = gsl_rng_uniform_pos( rngtype);
  std::size_t j {2} ;
 
  while( u > f[j]){ ++j ;}
  assert( j >= 2);
  return j ;
}


@*1 {\bf assign chromosomes to offspring}.
\label{sec:assignchromstooffspring}


assign chromosomes to potential offspring; there are $N$ parent pairs
so $2N$ diploid individuals and $4N$ gene copies; then parent pair $i$
will have chromosomes labelled $4i, 4i+1, 4i + 2, 4i + 3$ to
distribute among its offspring; the ancestry will record the ancestral
relations of the gene copies, or chromosomes


@<assign chromosomes to offspring@>=@#
static void assignchromstooffspring( const std::size_t& i,  const std::size_t &x, std::vector< std::pair<std::size_t,  std::size_t> >& y )
{
  /* \newline  |i| is the parent  pair index from 0 to $N-1$;  each parent pair has 4 chromosomes; 
    need |y| as vector of pairs since will eventually remove some of them; 
    each potential offspring is pair of gene copies; since will eventually remove some of them; 
    $0 \le 4i \le 4N-4$  */
  for(std::size_t k = 0 ; k < x ; ++k){
    y.push_back( std::make_pair( (4*i) + (gsl_rng_uniform_pos(rngtype) < 0.5 ? 0 : 1),  (4*i) + (gsl_rng_uniform_pos(rngtype) < 0.5 ? 2 : 3) ) ) ; }
}


@*1 {\bf add generation to population ancestry}.
\label{sec:addgeneration}

add a new set of surviving offspring to population ancestry; the
ancestry records the ancestral relations of gene copies 

@<add a new generation@>=@#
static void addgeneration( std::vector< std::size_t >& tree, const std::vector<double>& vcmf )
{
  /* \newline $N$ parent pairs produce offspring; $2N$ diploid individuals in the population */
  std::size_t x {} ;
  std::size_t s {} ;
  std::vector< std::pair<std::size_t, std::size_t>> y {} ;
  y.clear() ;
  /* \newline   |CONST_N| parent pairs produce potential offspring  \S~\ref{sec:constants} */ 
  for( std::size_t i = 0 ; i < CONST_N; ++i){
    /* \newline |x| is number of diploid potential offspring produced by parent pair |i|; \S~\ref{sec:randomX} */
    x = randomX(vcmf);
    s += x ;
    y.reserve( y.size() + x);
    /* \newline \S~\ref{sec:assignchromstooffspring} */
    assignchromstooffspring(i, x, y); }

/* \newline |y| is the record of diploid potential offspring */ 
  assert( y.size() >= CONST_POP_SIZE);

/* \newline we sample surviving offspring by shuffling an index vector |V|
and recording the offspring corresponding to the first $N$ elements of the shuffled vector */


std::vector<std::size_t> V (s,0);
  std::iota( V.begin(), V.end(), 0);
  /* \newline  |rng|    \S~\ref{sec:rngs} */
  std::shuffle( V.begin(), V.end(), rng);
  /* \newline we want to add the ancestral relations of  $4N$ gene copies to the ancestry */
  tree.reserve( tree.size() + (4*CONST_N) );
   /* \newline  add new generation to tree */

/* \newline |CONST_POP_SIZE| is $2N$ the number of diploid individuals */
  for( std::size_t z = 0 ; z < CONST_POP_SIZE; ++z){
    tree.push_back( y[V[z]].first );
    tree.push_back( y[V[z]].second );}
}


@*1 {\bf a random sample}.
\label{sec:randomsample}


get a random sample of diploid individuals by sampling a random set of
levels 

@<a random sample@>=@#
static void randomsample( std::vector< std::size_t>& V,  std::vector< std::size_t  >& sample )
{
  /* \newline  pair is (size of block, level of  immediate ancestor of block) */

  /* \newline  sample levels ranges from 0 to $4N - 2$  */
  /* \newline  the elements of |V| are in  $\{0, 2, 4, \ldots, 4N - 2\}$  */

  std::shuffle( V.begin(), V.end(), rng);
  
  assert( sample.size() == CONST_SAMPLE_SIZE);
  std::copy( V.begin(), V.begin() + CONST_SAMPLE_SIZE, sample.begin() );
}


@*1 {\bf get level of immediate ancestor}.
\label{sec:getagi}

look up the immediate ancestor of the individual living on level
$\ell$ at time $g$


@<immediate ancestor@>=@#
static std::size_t getagi( const std::size_t &g, const std::size_t& l, const std::vector< std::size_t>& t )
{
  /* \newline $0 \le  i < 4N$ */
  /* \newline  each tree 'row' or generation is  with indexes $0, 1, 2, \ldots, 4N - 1$ */
  /* \newline |CONST_POP_SIZE| is $2N$ the number of diploid individuals */
  assert( ((2*g*CONST_POP_SIZE) + l ) < t.size() );
  return t[ (2*g*CONST_POP_SIZE) + l ] ;
}




@*1 {\bf check tree completeness}.
\label{sec:treeincomplete}



checking if sample tree is incomplete; the tree is complete if the
leaves have a common ancestor 

@<incomplete@>=@#
static bool treeincomplete(  std::vector< std::size_t  >& sample,   const std::vector<std::size_t >& ancestry)
{
  std::unordered_set< std::size_t> s {};
  s.clear();
  assert( s.empty());

/* \newline |sample| has elements from $\{0,2,4,\ldots, 4N - 2\}$ since there are $2N$ diploid individuals; the
indexes of diploid individuals; diploid individual on level  $\ell$ has
chromosomes on levels  $\ell$ and $\ell + 1$  */
 
  std::for_each( sample.begin(), sample.end(), [&s](const auto &y){ s.insert( y); s.insert(y+1); });
  std::size_t  g = (ancestry.size()/(2*CONST_POP_SIZE)) - 1 ;
  std::vector<std::size_t> a (2*CONST_SAMPLE_SIZE );
  assert( a.size() == s.size());

   while( (s.size() > 1) && (g > 0) ){
     a.resize( s.size());
     std::transform( s.begin(), s.end(), a.begin(), [&g, &ancestry](const auto &i){ return getagi(g, i, ancestry);});
     s.clear();
     assert( s.empty());
     /* \newline |s| being an unordered set only records unique elements */
     std::for_each( a.begin(), a.end(), [&s]( const auto &b){ s.insert(b);});
     assert( s.size() > 0);
     --g ; }

   return s.size() > 1 ;
}


@*1 {\bf get a complete sample tree}.
\label{sec:getcompletesampletree}


get a complete sample tree by adding to the ancestry and drawing new
samples until a complete sample tree is found



@<one complete sample tree@>=@#
static void getcompletesampletree( std::vector<std::size_t>& ancestry,  std::vector<std::size_t>& sampledlevels, std::vector<std::size_t>& __V,  const std::vector<double>& vcmf )
{
@#

  /* \newline  \textcolor{violet}{\small \sf  get a complete sample tree 
   \newline  sample levels until a complete tree is found 
   \newline |tree| is population tree} */

@#

  /* \newline \textcolor{violet}{\small \sf   start a new population from scratch}  */
  ancestry.clear();
  ancestry.reserve( 2*CONST_POP_SIZE);
  std::vector<std::size_t>( 2*CONST_POP_SIZE).swap(ancestry);
  /* \newline \textcolor{violet}{\small \sf  initialise $A_{\ell}(0) = \ell$ for $\ell = 0,1, \ldots, N-1$} */
  std::iota( ancestry.begin(), ancestry.end(), 0);

 /* \newline \S~\ref{sec:addgeneration} */
  addgeneration(ancestry, vcmf);
  /* \newline \S~\ref{sec:randomsample} */
  randomsample( __V, sampledlevels );
  std::size_t number_generation = 1 ;
  /* \newline \S~\ref{sec:treeincomplete} */
  while( treeincomplete( number_generation, sampledlevels,  ancestry) ){
    addgeneration( ancestry, vcmf);
    randomsample(__V, sampledlevels);
    ++ number_generation ; }
}


@*1 {\bf update the site-frequency spectrum}.
\label{sec:updatesfs}

update the site-frequency spectrum (branch lengths) $\ell_{i}^{N}(n)$
given current block sizes


@<update the site-frequency spectrum@>=@#
static void updatesfs( std::vector<double>& __ells, std::unordered_map< std::size_t, std::size_t>& __m )
{

 @#

  std::for_each( __m.begin(), __m.end(), [&__ells](const auto &x){ assert(x.second > 0); assert( x.second < 2*CONST_SAMPLE_SIZE);  ++__ells[0]; ++__ells[x.second]; });
}


@*1 {\bf update estimates $\widehat r_{i}^{N}(n)$ of $\EE{ \widetilde R_{i}^{N}(n)}$}.
\label{sec:updaterrs}


given branch lengths $\ell_{i}^{N}(n)$ update estimates $\widehat r_{i}^{N}(n)$
of $\EE{ \widetilde R_{i}^{N}(n)}$
\begin{displaymath}
\widehat r_{i}^{N}(n) \leftarrow \widehat r_{i}^{N}(n)  +  \frac{\ell_{i}^{N}(n) }{ \sum_{j=1}^{2n-1} \ell_{j}^{N}(n)  }
\end{displaymath}

@<update  $\widehat r_{i}^{N}(n)$@>=@#
static void updaterrs( std::vector<double>& __errs, const std::vector<double>& __ells  )
{
 @#

  assert( __ells[0] > 0);
  @#
  
  const double d = __ells[0] ;

@#
  std::transform( __ells.begin(), __ells.end(), __errs.begin(), __errs.begin(), [&d]( const auto &x, const auto &y){ return y +  (x/d);});
}


@*1 {\bf read a complete tree}.
\label{sec:readcompletetree}

given the leaves of a complete tree read the branch lengths off the
tree


@<read a complete tree@>=@#
static void readcompletetree( std::vector<std::size_t> & sl, std::vector<double>& errs,  const std::vector<std::size_t>& __ancestry )
{ 

@#

  std::unordered_map< std::size_t, std::size_t> m {};
  m.clear(); @#
  std::unordered_map< std::size_t, std::size_t> tmp {};
  tmp.clear() ; @#
  std::vector<double> ells (2*CONST_SAMPLE_SIZE);

 @#

  std::size_t __g = (__ancestry.size()/(2*CONST_POP_SIZE)) - 1;

/* \newline \textcolor{violet}{\small initialise the map |m| with $2n$
blocks of size 1 each labelled with the sampled levels } */

std::for_each( sl.begin(), sl.end(), [&m](const auto &x){ m[x] = 1; m[x + 1] = 1; });

@#

assert(m.size() == 2*sl.size());

  while( m.size() > 1){
  /* \newline \textcolor{violet}{\small \sf  |updatesfs|   \S~\ref{sec:updatesfs}} */
    updatesfs( ells, m);
    tmp.clear();
    assert(tmp.size() < 1);
    /* \newline \textcolor{violet}{ \small  read new (merged and continuing) blocks into |tmp|; use |getagi| \S~\ref{sec:getagi}}  */
    std::for_each( m.begin(), m.end(), [&tmp, &__g, &__ancestry]( const auto &x ){ tmp[ getagi(__g, x.first, __ancestry)] += x.second; });
    m.clear();
    assert(m.size() < 1);
    m=tmp;
    assert( m.size() == tmp.size() );
    --__g ; }

 /* \newline \textcolor{violet}{ \small \sf  |updaterrs| \S~\ref{sec:updaterrs} } */
  updaterrs( errs, ells);
}


@*1 {\bf approximate $\EE{ \widetilde R_{i}^{N}(n)}$}. 
\label{sec:estimate}


estimate $\EE{ \widetilde R_{i}^{N}(n)}$ for a given
number of |CONST_EXPERIMENTS| \S~\ref{sec:constants}

@<go ahead --- get $\EE{\widetilde R_{i}^{N}(n)}$@>=@#
static void estimate( )
{ @#

  std::vector<std::size_t> ancestry__ {} ;
  std::vector<std::size_t> sl__ (CONST_SAMPLE_SIZE) ;
  std::vector<std::size_t> V__ (CONST_POP_SIZE);
  std::iota( V__.begin(), V__.end(), 0);
  std::transform( V__.begin(), V__.end(), V__.begin(), [](const auto &x){ return (x*2); });

@#
  std::vector<double> vcmfalpha__ (CONST_CUTOFF + 1);
  std::vector<double> vcmfkappa__ (CONST_CUTOFF + 1);
  /* \newline |generatecmf|  \S~\ref{sec:generatecmf} */
  generatecmf( vcmfalpha__, vcmfkappa__ );

@#

  std::vector<double> errs__ (2*CONST_SAMPLE_SIZE); 
  int r = CONST_EXPERIMENTS + 1;
  while( --r > 0){
  /* \newline \S~\ref{sec:getcompletesampletree} */
    getcompletesampletree( ancestry__, sl__, V__,  vcmfalpha__ );
    /* \newline \S~\ref{sec:readcompletetree} */
    readcompletetree( sl__, errs__, ancestry__); }
    /* \newline return the estimates $\widehat r_{i}^{N}(n)$ summed over the experiments */
  std::for_each( errs__.begin(), errs__.end(), [](const auto &x){ std::cout <<x << '\n' ; });
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
/* \newline \S~\ref{sec:generatecmf} */
@<generate cmf@>@#
/* \newline \S~\ref{sec:randomX} */
@<get one random number of potential offspring@>@#
/* \newline \S~\ref{sec:assignchromstooffspring} */
@<assign chromosomes to offspring@>@#
/* \newline \S~\ref{sec:addgeneration} */
@<add a new generation@>@#
/* \newline \S~\ref{sec:randomsample} */
@<a random sample@>@#
/* \newline \S~\ref{sec:treeincomplete} */
@<incomplete@>@#
/* \newline \S~\ref{sec:getagi} */
@<immediate ancestor@>@#
/* \newline \S~\ref{sec:getcompletesampletree} */
@<one complete sample tree@>@#
/* \newline \S~\ref{sec:updatesfs} */
@<update the site-frequency spectrum@>@#
/* \newline \S~\ref{sec:updaterrs} */
@<update  $\widehat r_{i}^{N}(n)$@>@#
/* \newline \S~\ref{sec:readcompletetree} */
@<read a complete tree@>@#
/* \newline \S~\ref{sec:estimate} */
@<go ahead --- get $\EE{\widetilde R_{i}^{N}(n)}$@>@#

int main()
{

/* \newline |setup_rng|   \S~\ref{sec:rngs} */
  setup_rng( static_cast<std::size_t>( atoi(argv[1])));
  /* \newline |estimate|  \S~\ref{sec:estimate} */
  estimate(); 

  gsl_rng_free( rngtype);

return GSL_SUCCESS ;
}


@* {\bf conclusions and bibliography}.
\label{sec:concl}


We estimate
$\widehat \varrho_{i}^{N}(n) = \EE{ \EE{R_{i}^{N}(n) :
\mathcal{A}^{(N,n)}}}$ when the sample comes from a finite diploid
panmictic population of constant size; the population may be evolving
acccording to random recruitment success.  To estimate
$\widehat \varrho_{i}^{N}(n)$ we evolve the population forward in time
and record the ancestry of the entire population; when a random sample
with a complete tree is found the ancestry of the sample is fixed and
we can read the branch lengths off the fixed complete tree; this
process is repeated a given number of times and so an estimate of
$\widehat \varrho_{i}^{N}(n)$ is obtained.  The idea behind
$\widehat \varrho_{i}^{N}(n)$ is the simple fact that the gene copies
of a sample share an ancestry, a gene genealogy, even if the genealogy
is unknown. One would want to compare $\widehat \varrho_{i}^{N}(n) $
to $\EE{R_{i}^{N}(n)} $ (mean relative branch lengths evaluated
without conditioning on the ancestry) to see the effect on gene
genealogies of conditioning on complete sample trees.








\begin{thebibliography}{1}

\bibitem{Birkner2018}
{\sc Birkner, M., Liu, H., and Sturm, A.}
\newblock Coalescent results for diploid exchangeable population models.
\newblock {\em Electronic Journal of Probability 23}, none (Jan. 2018).

\bibitem{kernighan1988c}
{\sc Kernighan, B.~W., and Ritchie, D.~M.}
\newblock The {C} programming language, 1988.

\bibitem{knuth1994cweb}
{\sc Knuth, D.~E., and Levy, S.}
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.

\bibitem{schweinsberg03}
{\sc Schweinsberg, J.}
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl 106\/} (2003), 107--139.

\bibitem{tange22gnu_paral}
{\sc Tange, O.}
\newblock {\em GNU {Parallel} 20221122}, 2022.
\newblock Zenodo.

\end{thebibliography}



@
\end{document}
