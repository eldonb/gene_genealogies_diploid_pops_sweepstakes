%$ NAFN=xindbetacoal_using_integer_partitions
%$ echo 'const double dalpha = 1.0;' > $NAFN.hpp
%$ echo 'const long double ALPHA = static_cast<long double>(dalpha);' >> $NAFN.hpp
%$ echo 'const long double KAPPA = 2.0L;' >> $NAFN.hpp
%$ echo 'const long double C_C = 1000.0L;' >> $NAFN.hpp
%$ echo 'const long double MM = (2.0L + ((1.0L  +  powl(2.0L , 1.0L - KAPPA))/(KAPPA - 1.0L)))/2.0L ;'>>$NAFN.hpp
%$ echo 'const long double CA = ((KAPPA + 2.0L) +(KAPPA*KAPPA))/2.0L;' >>$NAFN.hpp
%$ echo 'const long double CB = powl(2.0L , KAPPA)*((KAPPA - 2.0L)*(KAPPA - 1.0L));'>>$NAFN.hpp
%$ echo 'const long double CKAPPA = 2.0L*(KAPPA > 2.0L ? CA/CB : 1.0L)/(MM*MM);'>>$NAFN.hpp
%$ echo 'const double dgamma = 1.0 ;' >> $NAFN.hpp
%$ echo 'const long double GAMMA = static_cast<long double>(dgamma) ;'>>$NAFN.hpp
%$ echo 'const long double BETA = static_cast<long double>(gsl_sf_beta(2. - dalpha, dalpha)*(dgamma < 1. ? gsl_sf_beta_inc(2.-dalpha, dalpha, dgamma) : 1.));'>>$NAFN.hpp
%$ echo 'const long double CKAG = (CKAPPA + ((ALPHA*C_C)*BETA/powl(MM,ALPHA)))/4.0L;'>>$NAFN.hpp
%$ echo 'const unsigned int SAMPLE_SIZE = 5 ;'>>$NAFN.hpp
%$ echo 'const int EXPERIMENTS = 1e2 ;' >> $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ ./a.out $(shuf -i 43433-3838382 -n1)$ > myouts
%$ sed '1d' myouts | awk '{M=100.0; print log($1/M) - log(1 - ($1/M))}'> logitresout
%$ seq 9 | awk '{S=10;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ cweave $NAFN.w
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ emacs --version | head -n1 > innleggemacs
%$ g++ --version | head -n1 > innleggcpp
%$ lualatex --version | head -n2 > innlegglualatex
%$ cweave  --version | head -n1 > innleggcweave
%$ ctangle  --version | head -n1  > innleggctangle
%$ uname  --kernel-release -o  > innleggop
%$ bash --version | head -n1 > innleggbash
%$ sed -i 's/x86/x86\\/g' innleggbash
%$ gsl-config --version > innlegggsl
%$ spix --version > innleggspix
%$ parallel  --version | head -n1 > innleggparallel
%$ lualatex $NAFN.tex
%$ lualatex $NAFN.tex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright (C) 2025 aber der voll idiot
\documentclass[a4paper,10pt,final]{cweb}
\usepackage{fancyvrb,  graphicx, url}
\usepackage[svgnames]{xcolor}
\usepackage[inputlexerlinenos=true]{minted}
\setminted[cpp]{breaklines}
\usepackage{piton}
%% width=min
%% width=12cm
\PitonOptions{width=min,line-numbers,break-lines,indent-broken-lines,background-color=gray!5}
%% a4 paper size 210 x 297 millimeters
%% for xelatex
%%\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
%% \usepackage[all]{xy}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage{natbib}
\usepackage{abstract}
\renewcommand{\abstractname}{}
%\usepackage[backend=biber,style=authoryear-icomp,
%    sortlocale=de_DE,
%    natbib=true,
%    url=false, 
%    doi=true,
%    eprint=false]{biblatex}
%    \addbibresource{refs.bib}
    \usepackage[]{hyperref}
    %\hypersetup{colorlinks=true}
    \hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\usepackage[en-GB,showdow]{datetime2}
\DTMlangsetup[en-GB]{ord=raise,monthyearsep={,\space}}
\usepackage[rightcaption]{sidecap}
\sidecaptionvpos{figure}{m}
\usepackage{siunitx}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=newest}
\usepgfplotslibrary{units}
\sisetup{
  round-mode          = places,
  round-precision     = 3,
}
\setstretch{1}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\IN}{\ensuremath{\mathds{N}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\textbf{\small #1}}}
\newcommand{\s}{{\textcolor{blue}{\ensuremath\bullet}}}
\newcommand{\sa}[1]{\mathrel{\reflectbox{\rotatebox[origin=c]{#1}{\scalebox{0.5}{$\uparrow$}}}}}
\newcommand{\x}{{\circ}}
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
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  {\textsf{\textbf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\@@author}
\end{center}
\egroup
}
\makeatother
\title{Gene genealogies in diploid populations evolving according to
sweepstakes reproduction \\ --- approximating $\EE{R_i(n)}$ for the $\Omega$-$\delta_0$-Beta$(\gamma,2-\alpha,\alpha)$ coalescent}
 \author{Bjarki Eldon\footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}}  \footnote{ compiled @@
{\DTMcurrenttime} on  {\today}  \\
\input{innleggctangle} \\
 \input{innleggcpp} \\ kernel  \input{innleggop} \\  \input{innleggbash} \\
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  
\input{innlegglualatex}  \\ \input{innleggspix} \\  \input{innleggparallel} \\  \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391}}
\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

\begin{abstract} Let $\{ \xi^{n}(t) : t \ge 0 \}$ denote a coalescent,
$\# A$ the number of elements in a given finite set $A$, $n$ the
number of diploid individuals sampled (so that $2n$ gene copies or
chromosomes are sampled),
$L_{i}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in \xi^{n}(t) :
\#\xi = i \right\}dt $ and
$L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $ and
$\tau(n) \equiv \inf \left\{ t \ge 0 : \# \xi^{n}(t) = 1 \right\} $
for $i \in \{1, 2, \ldots, 2n-1\}$;
$R_{i}(n) \equiv L_{i}(n)/\sum_{j}L_{j}(n) $ for $i=1,2,\ldots, n-1$.
Then $L_{i}(n)$ is interpreted as the random total length of branches
supporting $i \in \{1, 2, \ldots, 2n-1\}$ leaves, with the length
measured in coalescent time units.  We then have
$L(n) = L_{1}(n) + \cdots + L_{2n-1}(n)$.  With this C++ code one
estimates the functionals $\EE{R_{i}(n)}$ when $\{\xi^n(t)\}$ is the
$\Omega$-$\delta_0$-Beta$(\gamma,2-\alpha,\alpha)$ coalescent where
$0 < \gamma \le 1$ and $0 < \alpha < 2$. The
$\Omega$-$\delta_0$-Beta$(\gamma,2-\alpha,\alpha)$ coalescents are a
family of $\Xi$-coalescents where the asynchronous merger sizes are
driven by a Beta$(\gamma,2-\alpha,\alpha)$ measure with an atom at 0;
the blocks participating in a merger are split into four groups
(parent chromosomes) uniformly at random (with replacement) and the
blocks in the same group are merged.
\end{abstract}


\tableofcontents



@* {\bf Copyright}. 
\label{sec:copy}



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
the GNU General Public License for more details.  You should have
received a copy of the GNU General Public License along with this file
(see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.



@* {\bf compilation and output}. 
\label{sec:compy}


This CWEB \citep{knuth1994cweb} document (the {\tt .w} file) can be
compiled with {\tt cweave} to generate a {\tt .tex} file, and with
{\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c} C++ code
file. 

Use the shell tool {\tt spix} on the  script appearing before  the
preamble (the lines starting with  \%\$); simply 

{\tt spix /path/to/the/sourcefile}

where {\tt sourcefile} is the {\tt .w} file

One   may also  copy the script into a file and run  {\tt parallel}
\cite{tange11:_gnu_paral} :

{\tt parallel --gnu -j1 :::: /path/to/scriptfile}



@* {\bf intro}.
\label{sec:intro}


Suppose a population is and evolves as in Definition
\begin{defn}
\label{sec:evol}
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



\begin{remark}
\label{sec:ill}
The mechanism described in Definition~\ref{sec:evol} is illustrated
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




Write $R_{i}(n) \equiv L_{i}(n)/\sum_{j=1}^{2n-1}L_{j}(n)$.  We
approximate $\EE{ R_{i}^{N}(n) }$ when the sample comes from a finite
diploid panmictic population evolving as in Definition~\ref{sec:evol}.


Let $X_{1}, \ldots, X_{N}$ denote the random number of potential
offspring produced by the $N$  current parent pairs. The
$X_{1},\ldots, X_{N}$  are independent.   Let $X$ be independent
of $X_{1}, \ldots, X_{N}$ and suppose
\begin{equation}
\label{eq:1}
\prb{X=k} =  C\left( \frac{1}{k^{a}} -  \frac{1}{(1+k)^{a}}  \right), \quad k \in \set{2,3, \ldots, \zeta(N)}
\end{equation}
with $a > 0$ and we choose $C$ such that  $\prb{X \in \{2,3, \ldots, \zeta(N) \} } = 1$.
Write $X\vartriangleright L(a,\zeta(N))$ when $X$ is distributed
according to \eqref{eq:1} for some given $a$ and $\zeta(N)$.  In any
given generation the current individuals independently produce
potential offspring according to \eqref{eq:1}; from the pool of
potential offspring $N$ of them are sampled uniformly at random and
without replacement to survive and reach maturity and replace the
current individuals. Note that \eqref{eq:1} is an extension of
\cite[Equation~11]{schweinsberg03}. Suppose  $X_{1},\ldots, X_{N}$ are
iid copies of  $X$ where   $\prb{X\vartriangleright
L(\alpha,\zeta(N))} = \varepsilon_{N}$, $\prb{X\vartriangleright
L(\kappa,\zeta(N))} = 1 - \varepsilon_{N} $ where $1 \le  \alpha < 2$
and $\kappa \ge 2$ both fixed.  


$r = \sum_{i}\one{k_{i}\ge 2} $ and $k \sum_{i}\one{k_{i}\ge 2}k_{i} $
\begin{subequations}
\begin{align}
%%\label{eq:3}
\lambda_{m;k_{1},\ldots,k_{r};s} &  = \binom{m}{k_{1}\cdots k_{r}s} \frac{1}{\prod_{j=2}^{m}\left( \sum\nolimits_{i}\one{k_{i}=j} \right)! }  \lambda_{n;k_{1},\ldots, k_{r};s}^{\prime}  \label{eq:5}  \\
\lambda_{n;k_{1},\ldots, k_{r};s}^{\prime} & = \one{k=2}\frac{C_{\kappa}}{C_{\kappa,\alpha,\gamma}} + \frac{c\alpha   }{m_{\infty}^{\alpha} C_{\kappa,\alpha,\gamma} } \sum_{\ell=0}^{s \wedge (4-r)} \binom{s}{\ell} \frac{(4)_{r + \ell}}{ 4^{k+\ell} }B(\gamma, k + \ell - \alpha, m - k - \ell + \alpha) \label{eq:4} \\
m  &  =  \frac 12 \left(2 +  \frac{1 + 2^{1 - \kappa} }{\kappa - 1}  \right)  \\ 
\gamma & = \one{\frac{\zeta(N) }{N} \to K } + \one{ \frac{\zeta(N) }{N} \to \infty  } \\
C_{\kappa, \alpha, \gamma} & =  \frac 14 C_{\kappa} +  \frac{\alpha c }{4 m^{\alpha} }B(\gamma,2-\alpha, \alpha) \\
C_{\kappa} & = \one{\kappa = 2} \frac{2 }{m^{2}} +  \one{\kappa > 2} \frac{2}{m^{2}} \frac{c_{\kappa}}{2^{\kappa}(\kappa - 2)(\kappa - 1)}
\end{align}
\end{subequations}
where $\kappa + 2 < c_{\kappa} < \kappa^{2}$ when $\kappa > 2$. In
\eqref{eq:4} $0 < \alpha < 2$. 


 The code follows in \S~\ref{sec:includes} -- \S~\ref{sec:main}; we
conclude in \S~\ref{sec:concs}.  Comments within the code are in \aths{this font and
color}




@* {\bf code}.
\label{sec:code}





@*1 {\bf the includes}.
\label{sec:includes}

the included libraries

@<includes@>=@#
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <numeric>
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
#include <unordered_set>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <fenv.h>
#include <unistd.h>
#include <limits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "xindbetacoal_using_integer_partitions.hpp"


@*1 {\bf the random number generators}.
\label{sec:rngs}


define the random number generators


\begin{Piton}
  std::random_device randomseed;
  std::mt19937_64 rng(randomseed());

gsl_rng * rngtype ;
static void setup_rng(const unsigned long int s)
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype, s) ;
}
\end{Piton}


@<rngs@>=@#
/* \newline \aths{ obtain a seed for the  random number engine}  */
  std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine}  */
  std::mt19937_64 rng(randomseed());


gsl_rng * rngtype ;
static void setup_rng( const  unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


@*1 {\bf exponential function}.
\label{sec:veldi}

compute $e^{x}$ checking for flows

\begin{Piton}
long double veldi( const long double v )
{
  feclearexcept(FE_ALL_EXCEPT);
  const long double svar =  expl( v ) ;
  return( fetestexcept( FE_OVERFLOW ) != 0 ? LDBL_MAX : ( fetestexcept(FE_UNDERFLOW) != 0 ? 0.0L : svar) ) ;
  
}
\end{Piton}


@<$e^{x}$@>=@#
long double veldi( const long double v )
{
  feclearexcept(FE_ALL_EXCEPT);
  const long double svar =  expl( v ) ;
  return( fetestexcept( FE_OVERFLOW ) != 0 ? LDBL_MAX : ( fetestexcept(FE_UNDERFLOW) != 0 ? 0.0L : svar) ) ;
}


@*1 {\bf the logarithm of the beta function}.
\label{sec:betafunc}

the logarithm of the beta$(x,a,b) \equiv   \int_{0}^{x} t^{a-1}(1-t)^{b-1}dt$

 the GSL incomplete beta function is normalised by the complete beta
 function

   the standard way would be | gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) | 

   $\one{x < 1}\left( \log f  +a\log x + b\log(1-x) - \log a \right) +
   \one{x = 1} \left(\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b) \right)    $ 


\begin{Piton}
static long double betafunc( const double x, const double a, const double b )
{
  assert( x > 0.);
  assert( a > 0.);
  assert( b > 0.);

  const long  double f = static_cast<long double>( (x < 1. ? gsl_sf_hyperg_2F1(a + b, 1, a+1.,x) : 1.) );  
  assert( lessthan(0.0L, f) );

  return( x < 1. ?  (logl( f ) +  ( static_cast<long double>( (a*log(x)) + (b * log(1-x)) - log(a) ) ) ) : (lgammal( static_cast<long double>(a)) + lgammal( static_cast<long double>(b) ) - lgammal( static_cast<long double>(a + b) ) ) ); 
}
\end{Piton}


@<$\log$ Beta$(x,a,b)$@>=@#
static long double betafunc( const double x, const double a, const double b )
{
  /* \newline  \aths{ the GSL incomplete beta function is normalised by the complete beta function}  */
  /* \newline \aths{ $0 < x \le 1$  is the cutoff point}  */
  assert( x > 0.);
  assert( a > 0.);
  assert( b > 0.);
  /* \newline \aths{ the standard way would be | gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) |} */
  /* \newline \aths{ return  the logarithm  of the beta function  as $\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b)$} */
  const long  double f = static_cast<long double>( (x < 1. ? gsl_sf_hyperg_2F1(a + b, 1, a+1.,x) : 1.) );
  


  /* \newline  \aths{ return $\one{x < 1}( \log f  +a\log x +
  b\log(1-x) - \log a) + \one{x = 1} (\log\Gamma(a) + \log\Gamma(b) -
  \log\Gamma(a+b))$}  */
  return( x < 1. ?  (logl( f ) +  ( static_cast<long double>( (a*log(x)) + (b * log(1-x)) - log(a) ) ) ) : (lgammal( static_cast<long double>(a)) + lgammal( static_cast<long double>(b) ) - lgammal( static_cast<long double>(a + b) ) ) ); 
}



@*1 {\bf number of collisions}.
\label{sec:numbercollisions}


compute $\prod_{j=2}^{n}\left( \sum_{i}\one{k_{i}} = j  \right)! $
where $k_{i}$ are the number of blocks  assigned to each  of 4  groups
(uniformly at random with replacement), and $n$ is the current number
of blocks, so that $2 \le \sum_{i}k_{i} \le n$ 

\begin{Piton}
static long double numbercollisions( const std::vector<unsigned int>& __k )
{
  assert( std::all_of( __k.cbegin(), __k.cend(), [](const auto x){return x > 1;}) );

  double l {} ;
  switch( __k.size() ){
  case 1 : {
    l = 1;  break ; }
  case 2 : {
    assert( __k[0] <= __k[1]); 
    l = (__k[0] == __k[1] ? 2 : 1);
    break ; }
  case 3 : {
    assert( __k[1] <= __k[2]);
    assert( __k[0] <= __k[1]);
    l = ( __k[0] == __k[2] ? 6 : ( __k[0] == __k[1] ? 2 : ( __k[1] == __k[2] ? 2 : 1)));
    break ; }
  case 4 : {
    assert( __k[2] <= __k[3]);
    assert( __k[1] <= __k[2]);
    assert( __k[0] <= __k[1]);
    l = ( __k[0] ==  __k[3] ? 24 : ( __k[0] == __k[2] ? 6 :  ( ( __k[0] == __k[1] ? ( __k[2] == __k[3] ? 4 : 2) : ( __k[1] == __k[3] ? 6 : ( __k[1] == __k[2] ? 2 : ( __k[2] == __k[3] ? 2 :1)))))));
    break ; }
  default : break ; }

  assert(l > 0) ;

  return static_cast<long double>(l) ;
}
\end{Piton}


@<product of collision numbers@>=@#
static long double numbercollisions( const std::vector<unsigned int>& __k )
{

@#

  assert( std::all_of( __k.cbegin(), __k.cend(), [](const auto x){return x > 1;}) );

@#

  double l {} ;

@#

  switch( __k.size() ){
  case 1 : {
    l = 1;  break ; }
  case 2 : {
    assert( __k[0] <= __k[1]); 
    l = (__k[0] == __k[1] ? 2 : 1);
    break ; }
  case 3 : {
    assert( __k[1] <= __k[2]);
    assert( __k[0] <= __k[1]);
    l = ( __k[0] == __k[2] ? 6 : ( __k[0] == __k[1] ? 2 : ( __k[1] == __k[2] ? 2 : 1)));
    break ; }
  case 4 : {
    assert( __k[2] <= __k[3]);
    assert( __k[1] <= __k[2]);
    assert( __k[0] <= __k[1]);
    l = ( __k[0] ==  __k[3] ? 24 : ( __k[0] == __k[2] ? 6 :  ( ( __k[0] == __k[1] ? ( __k[2] == __k[3] ? 4 : 2) : ( __k[1] == __k[3] ? 6 : ( __k[1] == __k[2] ? 2 : ( __k[2] == __k[3] ? 2 :1)))))));
    break ; }
  default : break ; }

  assert(l > 0) ;

  return static_cast<long double>(l) ;
}



@*1 {\bf descending factorial}.
\label{sec:ff}

descending factorial $(x)_{m} \equiv x(x-1)\cdots (x-m+1)$ and
$(x)_{0} \equiv 1$; return $\log (x)_{m}$

\begin{Piton}
long double ff( const long double x, const long double m)
{
  long double f = 1.0L;
  long double j = 0.0L;
  while(j < m){
    f *= (x - j);
    ++j ;
  }
  assert(f > LDBL_EPSILON );
  return logl(f) ;
}
\end{Piton}


@<$(x)_{m}$@>=@#
long double ff( const long double x, const long double m)
{
  long double f = 1.0L;
  long double j = 0.0L;
  while(j < m){
    f *= (x - j);
    ++j ;
  }
  assert(f > LDBL_EPSILON);
  return logl(f) ;
}


@*1 {\bf multinomial coefficient}.
\label{sec:multinomialconstant}

compute
\begin{displaymath}
\log \binom{n}{k_{1}\cdots k_{r} s} - \log \prod_{j=2}^{n} \left( \sum\nolimits_{i}\one{k_{i} = j} \right)!
\end{displaymath}
\begin{Piton}
long double multinomialconstant( const unsigned int m, const std::vector<unsigned int>& v_k)
{
  const long double s = lgammal( static_cast<long double>(1 +  m - std::accumulate( v__k.cbegin(), v__k.cend(), 0)) ) ;
  long double d =   static_cast<long double>(std::accumulate( v__k.cbegin(), v__k.cend(), 0., [](long double a, const auto x){return a + lgammal(static_cast<long double>(x+1)); })) ;
  return ( lgammal(static_cast<long double>(m + 1)) - d - s - logl(numbercollisions(v_k)) ) ;
}
\end{Piton}


@<$\log \tbinom{n}{k_{1}\cdots k_{r} s}$@>=@#
long double multinomialconstant( const unsigned int m, const std::vector<unsigned int>& v__k)
{
@#
  const long double s = lgammal( static_cast<long double>(1 +  m - std::accumulate( v__k.cbegin(), v__k.cend(), 0)) ) ;
  @#
  long double d =   static_cast<long double>(std::accumulate( v__k.cbegin(), v__k.cend(), 0.0L, [](long double a, const auto x){return a + lgammal(static_cast<long double>(x+1)); })) ;

@#
/* \newline \aths{\S~\ref{sec:numbercollisions}} */
return ( lgammal(static_cast<long double>(m + 1)) - d - s - logl(numbercollisions(v__k)) ) ;
}


@*1 {\bf $\lambda_{n;k_{1},\ldots, k_{r};s}$}.
\label{sec:lambdanks}

compute  $\lambda_{n;k_{1},\ldots, k_{r};s}$    \eqref{eq:5}


@<$\lambda_{n;k_{1},\ldots,k_{r};s}$@>=@#
long double lambdanks( const unsigned int m, const std::vector<unsigned int>& v__k )
{

  const unsigned int  r = v__k.size();
  const unsigned int k = std::accumulate( v__k.cbegin(),  v__k.cend(),0);
  assert( k > 1) ;
  assert( k <= m);
  const unsigned int s = m - k ;

  auto logchoose = [](const unsigned int x, const unsigned int y){ return (lgammal( static_cast<long double>(x+1)) - lgammal(static_cast<long double>(y + 1)) - lgammal(static_cast<long double>(x-y + 1)));};

  long double l = 0.0L;
  
  for( std::size_t ell = 0; ell <=  (s < 4-r ? s : 4-r); ++ell ){
    l += veldi(logchoose(s,ell)  +  ff(4.0L, static_cast<long double>(r+ell)) + betafunc(GAMMA, static_cast<long double>(k + ell) - ALPHA, static_cast<long double>(m-k-ell)+ALPHA) - (logl(4.0L)*static_cast<long double>(k+ell)) ); }

  l *= (ALPHA * C_C) ;
  l /= (CKAG * powl( MM, ALPHA) );

  l += (k < 3 ?  CKAPPA / CKAG : 0.0L);
  assert( l > LDBL_EPSILON ) ;

/* \newline \aths{\S~\ref{sec:multinomialconstant}} */
   return ( veldi( logl(l) +  multinomialconstant( m, v__k ) ) ) ;
}



@*1 {\bf structure {\tt RM\_t}}.
\label{sec:RM_t}

for storing the integer partitions (ordered merger sizes)


@<structure rmt @>=@#
struct RM_t {

@#

  std::vector< std::vector< unsigned int > > r__p {} ; 
} ;



@*1 {\bf generate integer partitions}.
\label{sec:GenPartitions}

generate ordered integer partitions of at least 2 elements and summing
to   |myInt| 

@<size 2 or larger integer partitions@>=@#
void GenPartitions( const unsigned int m,  std::vector<long double>& v__lambdan, std::vector< RM_t >& rates_partitions, std::vector< std::vector< long double > > & rates_for_sorting,  const unsigned int myInt,
                   const unsigned int PartitionSize,
                   unsigned int MinVal,
                   unsigned int MaxVal)
{

  long double rate {} ;
    std::vector<unsigned int> partition(PartitionSize);
    unsigned int idx_Last = PartitionSize - 1;
    unsigned int idx_Dec = idx_Last;   
    unsigned int idx_Spill = 0;        
    unsigned int idx_SpillPrev;        

    unsigned int LeftRemain = myInt - MaxVal - (idx_Dec - 1)*MinVal;  
    partition[idx_Dec] = MaxVal + 1; 

    

    do {
        unsigned int val_Dec = partition[idx_Dec] - 1;   
        partition[idx_Dec] = val_Dec;                    

        idx_SpillPrev = idx_Spill;         
        idx_Spill = idx_Dec - 1;           

        while (LeftRemain > val_Dec)       
        {
            partition[idx_Spill--] = val_Dec;
            LeftRemain -= val_Dec - MinVal;
           
        }  

        partition[idx_Spill] = LeftRemain; 
        

        char a = (idx_Spill) ? ~((-3 >> (LeftRemain - MinVal)) << 2) : 11; 
        char b = (-3 >> (val_Dec - LeftRemain));

        switch (a & b) 
        {
            case 1:
            case 2:
            case 3: idx_Dec = idx_Spill;
                    LeftRemain = 1 + (idx_Spill - idx_Dec + 1)*MinVal; 
                    break;

            case 5: for (++idx_Dec, LeftRemain = (idx_Dec - idx_Spill)*val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= MinVal); idx_Dec++)
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;

            case 6:
            case 7:
            case 11:idx_Dec = idx_Spill + 1;
                    LeftRemain += 1 + (idx_Spill - idx_Dec + 1)*MinVal;
                    break;


            case 9: for (++idx_Dec, LeftRemain = idx_Dec * val_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec + 1)); idx_Dec++)  
                        LeftRemain += partition[idx_Dec];

                    LeftRemain += 1 - (idx_Dec - 1)*MinVal;
                    break;

            case 10:for (LeftRemain += idx_Spill * MinVal + (idx_Dec - idx_Spill)*val_Dec + 1, ++idx_Dec; (idx_Dec <= idx_Last) && (partition[idx_Dec] <= (val_Dec - 1)); idx_Dec++)   
                        LeftRemain += partition[idx_Dec];

                    LeftRemain -= (idx_Dec - 1)*MinVal;
                    break;
        }

        while (idx_Spill > idx_SpillPrev)   
            partition[--idx_Spill] = MinVal;   
	
	const  unsigned  int  sama_summa   =    static_cast<unsigned int>(std::accumulate( partition.cbegin(), partition.cend(), 0)) == myInt ;
	
	 switch(sama_summa){
	  case 1 : {
        /* \newline \aths{\S~\ref{sec:lambdanks}} */        
        rate = lambdanks(m, partition);
	v__lambdan[m] += rate;
	(rates_partitions[m]).r__p.push_back( partition ) ;
	rates_for_sorting[m].push_back(rate) ;
	break ; }
	  case 0 : break ;
	  default : break; }
    } while (idx_Dec <= idx_Last);
}



@*1 {\bf all merger sizes when $m$ blocks}.
\label{sec:allmergers_for_m_blocks}





@<mergers for $m$ blocks@>=@#
static void allmergers_for_m_blocks( const unsigned int m, std::vector<long double>& v__lambdan, std::vector< RM_t >& rates_partitions, std::vector< std::vector< long double > > & rates_for_sorting )
{

  std::vector< unsigned int > onemerger (1);
  long double lm {} ;
  for( unsigned int s = 2; s <= m; ++s)
    {
      onemerger[0] = s ;
      /* \newline \aths{\S~\ref{sec:lambdanks}} */
      lm = lambdanks(m, onemerger);
      v__lambdan[m] += lm ;
      rates_partitions[m].r__p.push_back(  onemerger) ;
      rates_for_sorting[m].push_back(lm);
    }

  if( m > 3){
  /* \newline \aths{get the partitions of size at least 2} */
    unsigned int number_groups {} ;
    /* \newline \aths{$k$ is the sum of the merger sizes, partitions } */
    for( unsigned int k = 4; k <= m; ++k){
      number_groups =  k/2 > 4 ? 4 : k/2 ;
      for( unsigned int s = 2; s <= number_groups; ++s){
      /* \newline \aths{\S~\ref{sec:GenPartitions}} */
        GenPartitions(m, v__lambdan, rates_partitions, rates_for_sorting,  k, s, 2, k - (2*(s-1)) ); }
  } }
  assert( rates_for_sorting[m].size() > LDBL_EPSILON) ;
}



@*1 {\bf order the rates}.
\label{sec:order_rates}


order the rates \eqref{eq:5} in descending order for  generating CMFs for sampling
mergers 

@<put the rates in descending order@>=@#
static void order_rates(  const std::vector< long double > & v__rates_for_sorting, std::vector< unsigned int > & v__indx )
{
@#
  assert( v__rates_for_sorting.size() > 0) ;

@#
  v__indx.clear();

@#

  v__indx.resize( v__rates_for_sorting.size() ) ;

@#

  std::iota( v__indx.begin(), v__indx.end(), 0 ) ;

@#

  std::stable_sort(v__indx.begin(), v__indx.end(), [&v__rates_for_sorting](const unsigned int x, const unsigned int y){return v__rates_for_sorting[x] > v__rates_for_sorting[y];});
}



@*1 {\bf compute a CMF}.
\label{sec:generate_cmf_n_blocks}

generate a cumulative mass function from sorted rates \eqref{eq:5} for
sampling mergers

@<cmf@>=@#
static void generate_cmf_n_blocks( const long double lambdan,   const std::vector<long double>& v__rates, const std::vector<unsigned int> & v__indx, std::vector<long double>& v__cmf )
{

@#

  v__cmf.clear();

@#

  v__cmf.resize( v__rates.size() );

@#

  v__cmf[0] = v__rates[v__indx[0]] / lambdan ;

@#

  for( unsigned int i = 1; i < v__indx.size(); ++i)
    {
      v__cmf[i] = v__cmf[ i-1] + ( v__rates[ v__indx[i] ] / lambdan ); 
    }

@#

  assert( fabsl(v__cmf.back() - 1.0L) < 1.0e-9L ) ;
}



@*1 {\bf get all merger sizes}. 
\label{sec:allmergers}

generate all ordered  mergers

@<all mergers@>=@#
static void allmergers(  std::vector<long double>& v__lambdan, std::vector< RM_t >& a__rates_partitions, std::vector< std::vector< long double > > & a__rates_for_sorting,  std::vector< std::vector< unsigned int > > & a__indx, std::vector< std::vector< long double > > & a__cmfs)
{

@#

  for( unsigned int n = 2; n <= 2*SAMPLE_SIZE; ++n){
  /* \newline \aths{\S~\ref{sec:allmergers_for_m_blocks}} */
    allmergers_for_m_blocks(n, v__lambdan, a__rates_partitions, a__rates_for_sorting);
    assert(a__rates_for_sorting[n].size() > 0);
    assert(v__lambdan[n] > LDBL_EPSILON);
    /* \newline \aths{\S~\ref{sec:order_rates}} */
    order_rates( a__rates_for_sorting[n], a__indx[n]);
    /* \newline \aths{\S~\ref{sec:generate_cmf_n_blocks}} */
    generate_cmf_n_blocks(v__lambdan[n], a__rates_for_sorting[n], a__indx[n], a__cmfs[n]);
  }
}



@*1 {\bf udpate the lengths $\ell_{i}$}.
\label{sec:update_lengths}

@<update realisations of $L_{i}(n)$@>=@#
void update_lengths(const long double l,  const std::vector<unsigned int>& v__tree,   std::vector<long double>& v__lengths)
{

@#

  assert( l > LDBL_EPSILON);

@#

  const long double t = -logl(static_cast<long double>(1. - gsl_rng_uniform(rngtype))) / l ;

@#

  for( unsigned int i = 0; i < v__lengths.size(); ++i)
    {
      v__lengths[0] += t ;
      v__lengths[v__tree[i]] += t ;
    }
}



@*1 {\bf sample merger}.
\label{sec:pick_merger}


sample mergers

@<pick merger@>=@#
unsigned int pick_merger( const std::vector<long double>& v__cmf )
{
@#
  unsigned int j = 0;

@#

  const long double u = static_cast<long double>(gsl_rng_uniform(rngtype));
  assert( u <= 1.0L );
@#

  while( u > v__cmf[j]){++j;}

  return j ;
}




@*1 {\bf given mergers update tree}.
\label{sec:update_tree}

given merger sizes merge blocks and update tree


@<merge blocks and update tree@>=@#
void update_tree( std::vector<unsigned int> & v__tree, const std::vector<unsigned int>& merger_sizes )
{

@#

  assert( merger_sizes.size() > 0) ;

@#

  assert( static_cast<unsigned int>(std::accumulate( merger_sizes.cbegin(), merger_sizes.cend(),0)) <= v__tree.size() );

@#

  std::shuffle( v__tree.begin(), v__tree.end(), rng );

@#

  std::vector<unsigned int > newblocks {};

@#

  newblocks.clear();

@#

  newblocks.reserve( merger_sizes.size() ) ;

@#

  assert( newblocks.size() < 1) ;

@#

  for( const auto &s : merger_sizes){
    newblocks.push_back( std::accumulate( v__tree.crbegin(), v__tree.crbegin() + s, 0));
    v__tree.resize( v__tree.size() - s) ; }

@#

  v__tree.insert( v__tree.end(), newblocks.cbegin(), newblocks.cend() );
}





@*1 {\bf one experiment}.
\label{sec:one_experiment}


one experiment 

@<one realisation@>=@#
void one_experiment( const std::vector<long double>& v__lambdan,  std::vector<long double>&  v__lengths, const  std::vector< std::vector< long double > > & a__cmfs, const std::vector< RM_t >& a__partitions, const std::vector< std::vector< unsigned int> > & a__indx, std::vector<long double>& v__ri )
{

@#

  std::vector< unsigned int> v__tree( 2*SAMPLE_SIZE, 1);


@#

  std::fill( v__lengths.begin(), v__lengths.end(), 0.0L) ;

@#

  unsigned int m {} ;
  while( v__tree.size() > 1)
    {
      /* \newline \aths{\S~\ref{sec:update_lengths}} */
      update_lengths( v__lambdan[v__tree.size()], v__tree, v__lengths );
      /* \newline \aths{\S~\ref{sec:pick_merger}} */
      m = pick_merger( a__cmfs[v__tree.size()]  ) ;
      /* \newline \aths{\S~\ref{sec:update_tree}} */
      update_tree( v__tree, a__partitions[v__tree.size()].r__p[ a__indx[v__tree.size()][m] ] ) ;
    }
  assert( v__tree.back() == (2*SAMPLE_SIZE));

@#

  assert( v__lengths[0] > LDBL_EPSILON) ;

@#

  const long double d = v__lengths[0];

@#

  std::transform( v__lengths.cbegin(), v__lengths.cend(), v__ri.begin(), v__ri.begin(), [&d](const auto x, const auto y){ return y + (x/d);});
}

@*1 {\bf approximate $\EE{R_{i}(n)}$}.
\label{sec:approximate}


@<go ahead -- approximate $\EE{R_{i}(n)}$@>=@#
void approximate()
{

@#

  std::vector<long double> v__lambdan ((2*SAMPLE_SIZE) + 1) ;

@#

/* \newline \aths{\S~\ref{sec:RM_t}} */
  std::vector< RM_t > a__parts ( (2*SAMPLE_SIZE) + 1,  RM_t {});

@#

  std::vector< std::vector< long double > > a_rates_sorting (  (2*SAMPLE_SIZE) + 1, std::vector<long double> {}) ;

@#


  std::vector< std::vector< long double > > a__cmfs ( (2*SAMPLE_SIZE)  + 1, std::vector< long double > {}) ;

@#

  std::vector< std::vector< unsigned int > > a__indx ( (2*SAMPLE_SIZE)  + 1, std::vector< unsigned int > {}) ;

@#

  std::vector<long double> v__lengths (2*SAMPLE_SIZE) ;

@#

  std::vector< long double> v__ri (2*SAMPLE_SIZE) ;

@#

/* \newline \aths{\S~\ref{sec:allmergers}} */
  allmergers(v__lambdan, a__parts, a_rates_sorting, a__indx, a__cmfs) ;

@#

  int r = EXPERIMENTS + 1;

@#

  while( --r > 0)
  /* \newline \aths{\S~\ref{sec:one_experiment}} */
    { one_experiment(v__lambdan, v__lengths, a__cmfs, a__parts, a__indx, v__ri) ; }

  for( const auto &z:v__ri)
    {std::cout << z << '\n';}
}






@*1 {\bf main}.
\label{sec:main}

the |main| module

@C


/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:veldi}} */
@<$e^{x}$@>@#
/* \newline \aths{\S~\ref{sec:betafunc}} */
@<$\log$ Beta$(x,a,b)$@>@#
/* \newline \aths{\S~\ref{sec:numbercollisions}} */
@<product of collision numbers@>@#
/* \newline \aths{\S~\ref{sec:ff}} */
@<$(x)_{m}$@>@#
/* \newline \aths{\S~\ref{sec:multinomialconstant}} */
@<$\log \tbinom{n}{k_{1}\cdots k_{r} s}$@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<$\lambda_{n;k_{1},\ldots,k_{r};s}$@>@#
/* \newline \aths{\S~\ref{sec:RM_t}} */
@<structure rmt @>@#
/* \newline \aths{\S~\ref{sec:GenPartitions}} */
@<size 2 or larger integer partitions@>@#
/* \newline \aths{\S~\ref{sec:allmergers_for_m_blocks}} */
@<mergers for $m$ blocks@>@#
/* \newline \aths{\S~\ref{sec:order_rates}} */
@<put the rates in descending order@>@#
/* \newline \aths{\S~\ref{sec:generate_cmf_n_blocks}} */
@<cmf@>@#
/* \newline \aths{\S~\ref{sec:allmergers}} */
@<all mergers@>@#
/* \newline \aths{\S~\ref{sec:update_lengths}} */
@<update realisations of $L_{i}(n)$@>@#
/* \newline \aths{\S~\ref{sec:pick_merger}} */
@<pick merger@>@#
/* \newline \aths{\S~\ref{sec:update_tree}} */
@<merge blocks and update tree@>@#
/* \newline \aths{\S~\ref{sec:one_experiment}} */
@<one realisation@>@#
/* \newline \aths{\S~\ref{sec:approximate}} */
@<go ahead -- approximate $\EE{R_{i}(n)}$@>@#

int main(int argc, const char * argv[])
{

/* \newline \aths{\S~\ref{sec:rngs}} */
 setup_rng( static_cast< std::size_t>( atoi(argv[1])) ) ;

/* \newline \aths{\S~\ref{sec:approximate}} */

   approximate() ;


@#
  
  gsl_rng_free( rngtype);

  return 0 ;
}


@* {\bf conclusion  and bibliography}.
\label{sec:concs}

we approximate the functionals   $\EE{R_{i}(n)}$ for $i = 1,2,\ldots ,
2n-1$  when  the coalescent
is the  $\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$
coalescent with transition  rates as in \eqref{eq:5}.
Figure~\ref{fig:graph} holds  an example 


\begin{SCfigure}[0.8][htb]
    \begin{tikzpicture}
      \begin{axis}[
        xlabel = $\log(i/n) - \log(1 - i/n)$,
        axis line style = {draw = none},
        tick style = {draw=none},
        xticklabels={draw=none},
        yticklabels={draw=none},
        legend pos=north east]
        \addplot table[col sep=comma] {forplottingfile1};
        %\addplot table[col sep=comma] {forplottingfile2};
        %\addplot table[col sep=comma] {forplottingfile3};
        \addlegendentry{$\gamma = 0.1$}
        %\addlegendentry{$\gamma = 0.5$}
        %\addlegendentry{$\gamma = 1$}
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{ \it  An example approximation of $\EE{R_{i}(n)}$
 for the  given  parameter values and graphed as logits against
 $\log(i/n) - \log(1 - i/n)$ for $i = 1,2,\ldots, 2n-1$ where $n$ is
 sample size}}
       \label{fig:graph}
       \end{SCfigure}







\bibliographystyle{plain}
\begin{thebibliography}{99}

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
