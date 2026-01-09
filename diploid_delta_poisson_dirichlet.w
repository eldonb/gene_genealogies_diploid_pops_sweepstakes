%$ NAFN=diploid_delta_poisson_dirichlet
%$ echo 'const unsigned int SAMPLE_SIZE = 50;' > $NAFN.hpp
%$ echo 'const double ALPHA = 0.01;' >> $NAFN.hpp
%$ echo 'const double KAPPA = 2;' >> $NAFN.hpp
%$ echo 'const double CEPS = 100;' >> $NAFN.hpp
%$ echo 'const double MINF = (2. + (1 + pow(2., 1-KAPPA)/(KAPPA -1)))/2. ;' >> $NAFN.hpp
%$ echo 'const double BKAPPA = pow(2,KAPPA)*(KAPPA - 2)*(KAPPA - 1);' >> $NAFN.hpp
%$ echo 'const double AKAPPA = ((KAPPA+2) + pow(KAPPA,2))/2.;' >> $NAFN.hpp
%$ echo 'const double CKAPPA = ((KAPPA > 2 ? AKAPPA / BKAPPA : 1.) *2.) / pow(MINF,2.) ;' >> $NAFN.hpp
%$ echo 'const int EXPERIMENTS = 1e5 ;' >> $NAFN.hpp
%$ ctangle $NAFN.w
%$ NAFN=diploid_delta_poisson_dirichlet
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 323383-1919101019 -n1) | sed '1d' | awk '{S=1e5;print log($1/S) - log(1 - ($1/S))}' > logitresout
%$ seq 49 | awk '{S=50;print log($1/S) - log(1 - ($1/S))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ sed -i 's/ALPHA = 0.01/ALPHA = 0.99/g' $NAFN.hpp
%$ sed -i 's/CEPS = 100/CEPS = 1/g' $NAFN.hpp
%$ ctangle $NAFN.w
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
%$ rm -f gg_*_.txt
%$ ./a.out $(shuf -i 323383-1919101019 -n1) | sed '1d' | awk '{S=1e5;print log($1/S) - log(1 - ($1/S))}' > logitresout
%$ paste -d',' nlogits logitresout > forplottingfile2
%$ cweave $NAFN.w
%$ tail -n4 $NAFN.tex > endi
%$ for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
%$ cat endi >> $NAFN.tex
%$ emacs --version | head -n1 > innleggemacs
%$ g++ --version | head -n1 > innleggcpp
%$ xelatex --version | head -n1  > innleggxelatex
%$ cweave  --version | head -n1 > innleggcweave
%$ ctangle  --version | head -n1  > innleggctangle
%$ uname  --kernel-release -o  > innleggop
%$ bash --version | head -n1 > innleggbash
%$ sed -i 's/x86/x86\\/g' innleggbash
%$ gsl-config --version > innlegggsl
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
\PitonOptions{width=min,line-numbers,break-lines,indent-broken-lines,background-color=gray!5}
%% a4 paper size 210 x 297 millimeters
%% for xelatex
%%\usepackage{xunicode}
\usepackage{fontspec}
\usepackage{xltxtra}
\usepackage{lineno}
%%\usepackage[all]{xy}
%%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage{natbib}
%\usepackage[backend=biber,style=authoryear-icomp,
%    sortlocale=de_DE,
%    natbib=true,
%    url=false, 
%    doi=true,
%    eprint=false]{biblatex}
%    \addbibresource{refs.bib}
    \usepackage[]{hyperref}
 %%   \hypersetup{colorlinks=true}
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
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  {\textsf{\large \textbf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\large \@@author}
\end{center}
\egroup
}
\makeatother
%%\pagecolor{DarkGray}
\title{Gene genealogies in diploid populations evolving according to
sweepstakes reproduction \\ --- approximating $\EE{R_i(n)}$ for the
 $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent}
 \author{Bjarki Eldon
 \footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{ compiled @@
{\DTMcurrenttime} on  {\today}  \\
\input{innleggctangle} \\
 \input{innleggcpp} \\ kernel  \input{innleggop} \\  \input{innleggbash} \\
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  {\LaTeX}
\input{innleggxelatex}  \\  written using  \input{innleggemacs} }
\orcidlink{0000-0001-9354-2391} }


\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


\begin{abstract}
 Let $\# A$ denote the  number of elements in a finite
set $A$.     For a given coalescent $\set{\xi^{n}}$  write
$L_{i}(n) \equiv  \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :
\#\xi = i }dt$ and  $L(n) \equiv  \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $
where $\tau(n) \equiv  \inf \set{t \ge 0 : \#\xi^{n}(t) = 1}$. Then
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Write $R_{i}(n) \equiv
L_{i}(n)/L(n)$ for $i = 1,2,\ldots, n-1$.   With this C++ simulation  code we
approximate   $\EE{R_{i}(n)}$ when      $\set{\xi^{n}} \equiv \set{\xi^{n}(t) : t \ge 0
} $ is  the   $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$
coalescent with $0 < \alpha < 1$,  i.e.\ the $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$
coalescent derived from a diploid panmictic population of constant
size evolving according to a specific model of sweepstakes
reproduction (heavy-tailed offspring number distribution) without
selfing (each diploid offspring has 2 diploid parents).     
\end{abstract}

\tableofcontents



@* {\bf Copyright}.
\label{sec:copyright}

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
\label{sec:compile}


This CWEB \citep{knuth1994cweb} document (the {\tt .w} file) can be
compiled with {\tt cweave} to generate a {\tt .tex} file, and with
{\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c} C++ code
file. 


Use the shell tool {\tt spix} on the  script appearing before  the
preamble (the lines starting with  \%\$); simply 

{\tt spix /path/to/the/sourcefile}

where {\tt sourcefile} is the {\tt .w} file

\begin{Verbatim}[frame=leftline,numbers=left,fontsize=\small,fontshape=it,formatcom=\color{teal}]
NAFN=diploid_delta_poisson_dirichlet
echo 'const unsigned int SAMPLE_SIZE = 50;' > $NAFN.hpp
echo 'const double ALPHA = 0.01;' >> $NAFN.hpp
echo 'const double KAPPA = 2;' >> $NAFN.hpp
echo 'const double CEPS = 100;' >> $NAFN.hpp
echo 'const double MINF = (2. + (1 + pow(2., 1-KAPPA)/(KAPPA -1)))/2. ;' >> $NAFN.hpp
echo 'const double BKAPPA = pow(2,KAPPA)*(KAPPA - 2)*(KAPPA - 1);' >> $NAFN.hpp
echo 'const double AKAPPA = ((KAPPA+2) + pow(KAPPA,2))/2.;' >> $NAFN.hpp
echo 'const double CKAPPA = ((KAPPA > 2 ? AKAPPA / BKAPPA : 1.) *2.) / pow(MINF,2.) ;' >> $NAFN.hpp
echo 'const int EXPERIMENTS = 1e5 ;' >> $NAFN.hpp
ctangle $NAFN.w
NAFN=diploid_delta_poisson_dirichlet
g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
rm -f gg_*_.txt
./a.out $(shuf -i 323383-1919101019 -n1) P sed '1d' P awk '{S=1e5;print log($1/S) - log(1 - ($1/S))}' > logitresout
seq 49 P awk '{S=50;print log($1/S) - log(1 - ($1/S))}' > nlogits
paste -d',' nlogits logitresout > forplottingfile1
sed -i 's/ALPHA = 0.01/ALPHA = 0.99/g' $NAFN.hpp
sed -i 's/CEPS = 100/CEPS = 1/g' $NAFN.hpp
ctangle $NAFN.w
g++ -std=c++26 -m64 -march=native -O3 -x c++ $NAFN.c -lm -lgsl -lgslcblas
rm -f gg_*_.txt
./a.out $(shuf -i 323383-1919101019 -n1) P sed '1d' P awk '{S=1e5;print log($1/S) - log(1 - ($1/S))}' > logitresout
paste -d',' nlogits logitresout > forplottingfile2
cweave $NAFN.w
tail -n4 $NAFN.tex > endi
for i in $(seq 5); do $(sed -i '$d' $NAFN.tex) ; done
cat endi >> $NAFN.tex
emacs --version P head -n1 > innleggemacs
g++ --version P head -n1 > innleggcpp
xelatex --version P head -n1  > innleggxelatex
cweave  --version P head -n1 > innleggcweave
ctangle  --version P head -n1  > innleggctangle
uname  --kernel-release -o  > innleggop
bash --version P head -n1 > innleggbash
sed -i 's/x86/x86\\/g' innleggbash
gsl-config --version > innlegggsl
xelatex $NAFN.tex
\end{Verbatim}
where {\tt P} is the system pipe operator.   Figure~\ref{fig:graph} records an example
of estimates of $\EE{R_{i}(n)}$. 


One   may also  copy the script into a file and run  {\tt parallel}
\cite{tange11:_gnu_paral} :

{\tt parallel --gnu -j1 :::: /path/to/scriptfile}




@* {\bf intro}.
\label{sec:intro}




Let $n, r, k_{1},\ldots, k_{r},\in \IN$, $n,k_{1},\ldots, k_{r} \ge
2$,  $\sum_{i}k_{i} \le n$, $s = n - \sum_{i} k_{i}$.   
The $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent
$\set{\xi^{n}}$   has
transition rate
\begin{equation}
\label{eq:1}
\begin{split}
\lambda_{n;k_{1},\ldots, k_{r};s} &  = \one{r=1,k_{1}=2} \binom{n}{2} \frac{C_{\kappa}}{C_{\kappa} + c(1-\alpha) }  \\
 & +   \binom{n }{k_{1}\ldots k_{r}\, s} \frac{1}{\prod_{j=2}^{n} \left( \sum_{i}\one{k_{i} = j} \right)!  }   \frac{c}{C_{\kappa} + c(1-\alpha)} p_{n;k_{1},\ldots,k_{r};s}
\end{split}
\end{equation}
where $0 < \alpha < 1$, $\kappa \ge 2$, $c \ge 0$ all fixed, and 
\begin{equation}
\label{eq:2}
\begin{split}
 p_{n;k_{1},\ldots,k_{r};s}  & =  \frac{\alpha^{r + s - 1} \Gamma(r+s)  }{\Gamma(n) } \prod_{i=1}^{r}(k_{i} - \alpha - 1)_{k_{i} - 1} \\
 C_{\kappa} & =   \one{\kappa = 2} \frac{2}{m_{\infty}^{2}} + \one{\kappa > 2} c_{\kappa}
\end{split}
\end{equation}
and $\kappa + 2 < c_{\kappa} < \kappa^{2} $ for $\kappa > 2$.    The
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent is an example of
a simultaneous multiple-merger coalescent. 
Simultaneous mergers  in    up to $\lfloor n/2 \rfloor $ groups for  
 $n$  current number of  blocks   may occur at such times.  



The $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent is
a diploid version of the $\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$
coalescent.  In the diploid version, each group of blocks in a
`partition' sampled from the $\delta_0$-Poisson-Dirichlet$(\alpha,0)$
coalescent is split into 4 boxes uniformly at random and the blocks in
the same box are merged.  Since the
$\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent can be
obtained from a model of a diploid panmictic population evolving
absent selfing and according to sweepstakes reproduction the boxes
represent the 4 parent chromosomes (2 from each parent) involved in a
large offspring number event.  With this C++ code we approximate
$\EE{R_{i}(n)}$ when the coalescent is the
$\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ 

We sample the  $\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$
coalescent by {\it (i)}  sampling a partition from the
$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent and {\it (ii)} splitting
the  partition sizes (blocks)  into boxes uniformly at random and with
replacement; we repeat steps {\it (i)} and {\it (ii)} until we have at
least one box with at least 2 blocks.  A merger then occurs and we
merge blocks and update the tree (the current block sizes).   


The algorithm is summarised in \S~\ref{sec:code}, the code follows in
\S~\ref{sec:includes} -- \S~\ref{sec:main}; we conclude in \S~\ref{sec:concl}.   Comments within the code are in
\aths{this font and colour} 


@* {\bf code}.
\label{sec:code}


we summarise the algorithm;
\begin{equation}
\label{eq:3}
\lambda_{m} = 
 \sum_{\substack{ 2 \le k_{1} \le \ldots \le k_{r} \le  m \\ k_{1} + \cdots + k_{r}   \le m  }}
\lambda_{m;k_{1},\ldots, k_{r}; s}
\end{equation}
\begin{enumerate}
\item generate all possible (simultaneous ordered)   mergers for $m=
2,3,\ldots, n$ blocks  \S~\ref{sec:allpartitionsrates}
\item for every merger(s) sizes   compute and record the transition
rate \eqref{eq:1}
\item record the total transition rate $\lambda_{m}$ \eqref{eq:3}   for
$m=2,3,\ldots, m$ for sampling partitions 
\item $(r_{1},\ldots, r_{n-1}) \leftarrow (0,\ldots, 0)$
\item for each of $M$ experiments : \S~\ref{sec:approximateeri}
\begin{enumerate}
\item $(\ell_{1}, \ldots, \ell_{n-1}) \leftarrow (0,\ldots, 0) $
\item $m \leftarrow n$
\item $(\xi_{1},\ldots, \xi_{n}) \leftarrow (1,\ldots, 1) $
\item {\bf while}   $m > 1$ : \S~\ref{sec:oneexperiment}
\begin{enumerate}
\item $t \leftarrow $ Exp($\lambda_{m})$ 
\item $\ell_{\xi} \leftarrow  \ell_{\xi} + t $ for $\xi = \xi_{1},
\ldots, \xi_{m}$ 
\item sample a partition  $k_{1}, \ldots, k_{r}$ \S~\ref{sec:samplemerger}  using the transition
rates \eqref{eq:1} and split partition into boxes \S~\ref{sec:splitgroups}  until at least one
box has at least 2 blocks   \S~\ref{sec:untilmerger}
\item given merger sizes merge blocks \S~\ref{sec:updatetree}
\item $m \leftarrow m - \sum_{i}b_{i}\one{b_{i} > 1}  +
\sum_{i}\one{b_{i} > 1} $ where $b_{i}$ is the number of blocks in box $i$
\end{enumerate}
\item $r_{i} \leftarrow  \ell_{i}/\sum_{j}\ell_{j} $ for  $i =
1,2,\ldots, n-1$
\end{enumerate}
\item return $r_{i}/M$ as an approximation of $\EE{R_{i}(n)}$ for $i=1,2,\ldots,n-1$
\end{enumerate}




@*1 {\bf includes}.
\label{sec:includes}

the included libraries

@<includes@>=@#
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <random>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <assert.h>
#include <float.h>
#include <fenv.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <boost/math/special_functions/factorials.hpp>
#include "diploid_delta_poisson_dirichlet.hpp"





@*1 {\bf random number generators}.
\label{sec:rngs}

define the STL and GSL random number generators 


the STL rngs 

%%\begin{minted}[bgcolor=LightGray,linenos,breaklines,mathescape]{cpp}
\begin{Piton}[background-color=gray!5,line-numbers,width=min,auto-gobble]
std::random_device randomseed;
  std::mt19937_64 rng(randomseed());
\end{Piton}
%%\end{minted}


the GSL random number generator 


%%\begin{minted}[bgcolor=LightGray,linenos,breaklines,mathescape]{cpp}
\begin{Piton}[background-color=gray!5,line-numbers,width=min,auto-gobble]
gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}
\end{Piton}
%%\end{minted}


@<rngs@>=@#
/* \newline \aths{ obtain a seed out of thin air for the random number engine} */
  std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine} */
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




@*1 {\bf descending factorial}.
\label{sec:descfact}

compute the descending factorial $(x)_{m} \equiv x(x-1)\cdots (x-m+1)$
and $(x)_{0} \equiv 1$.

using the boost library:
%% [background-color=gray!5,line-numbers,width=min,auto-gobble]
\begin{Piton}
static double descending_factorial(const double x, const double m)
{
 return static_cast<double>( boost::math::falling_factorial( x, static_cast<unsigned int>(m)) ) ;
 }
\end{Piton}

A direct implementation is
\begin{Piton}
 double p = 1;
  for( double i = 0; i < m; ++i){
    p *= (x - i); }
  return p ;
\end{Piton}



@<compute descending factorial@>=@#
static double descending_factorial(const double x, const double m)
{

 /* \newline \aths{using the boost library function : } */
 return static_cast<double>( boost::math::falling_factorial( x, static_cast<unsigned int>(m)) ) ;
  
/* \newline \aths{direct implementation : }
 | double p = 1;
  for( double i = 0; i < m; ++i){
    p *= (x - i); }
  return p ; |
  */
}



@*1 {\bf a power function}.
\label{sec:power}


compute $x^{y}$ checking for under and overflow 


\begin{Piton}
static double veldi( const double x, const double y )
  {
    feclearexcept(FE_ALL_EXCEPT);
    const double d = pow(x,y);

    return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;   
  }
\end{Piton}


@<checked power function@>=@#
static double veldi( const double x, const double y )
  {
    feclearexcept(FE_ALL_EXCEPT);
    const double d = pow(x,y);

    return( fetestexcept(FE_UNDERFLOW) ? 0. : (fetestexcept(FE_OVERFLOW) ? FLT_MAX : d) ) ;   
  }


@*1 {\bf $\lambda_{n;k_{1},\ldots, k_{r};s}$ \eqref{eq:1}}.
\label{sec:lambdanks}


compute the rate $\lambda_{n;k_{1},\ldots, k_{r};s}$ \eqref{eq:1}  of
the (ordered)   partition $(k_{1},\ldots,k_{r})$ 


@<compute  $\lambda_{n;k_{1},\ldots, k_{r};s}$ \eqref{eq:1}@>=@#
static double lambdanks( const double n, const std::vector<unsigned int>& v_k )
{

@#

/* \newline \aths{|v_k| is the partition to be split into boxes} */

  assert( v_k.size() > 0);

@#

  assert(std::all_of( v_k.cbegin(), v_k.cend(), [](const auto k){return k > 1;}));

@#
  
  double d {};
  double k {} ;
  double f {1} ;
  const double r = static_cast<double>(v_k.size());

@#

 /* \newline \aths{the counts $\sum_{i}\one{k_{i} = j}$ for
 $j=2,\ldots, n$  } */
  std::unordered_map<unsigned int, unsigned int> counts {} ;
  
  for( std::size_t i = 0; i < v_k.size(); ++i){
    /* \newline  \aths{ $f = \prod_{i=1}^{r} (k_{i} - \alpha - 1)_{k_{i} - 1}$} */
    f *= descending_factorial( static_cast<double>(v_k[i]) - 1. - ALPHA,  static_cast<double>(v_k[i]) -1);
    /* \newline  \aths{count number of  occurrences $c_{j} = \sum_{i}\one{k_{i}=j} $  of each block in
    the partition $(k_{1},\ldots , k_{r})$ } */
    ++counts[v_k[i]];
    /* \newline \aths{$k =  k_{1} + \cdots + k_{r}$ } */
    k += static_cast<double>(v_k[i]) ;
    /* \newline \aths{$d = \sum_{i} \log\Gamma(k_{i}+1) $ } */
    d += lgamma(static_cast<double>(v_k[i] + 1)) ; }

@#

  assert( k < n + 1 );

@#

  const double s = n - k;

@#

  /* \newline \aths{ $p = \sum_{j} \log\Gamma(c_{j} + 1)$ } */
  const double p = static_cast<double>(std::accumulate( counts.begin(), counts.end(), 0, [](double a, const auto &x){ return a +  lgamma( (double)x.second + 1);}));

@#

 
  const double l = ((v_k.size() < 2 ? (v_k[0] < 3 ? 1. : 0) : 0)*CKAPPA) + (CEPS * veldi(ALPHA, r+s-1) * tgamma(r+s)   *f/tgamma(n) ) ;


@#

  return ( veldi(exp(1), (lgamma( n + 1.) - d)  - lgamma( n - k + 1 ) - p) * l / (CKAPPA + (CEPS*(1-ALPHA))) ) ;
}



@*1 {\bf generate all fixed-size partitions  summing to a given number}.
\label{sec:partsumm}

generate all  partitions  of fixed size  |PartitionSize|   summing to $|myInt| \le m $ and compute the rate
$\lambda_{m;k_{1},\ldots, k_{r};s}$ for it when $m$ is the number of blocks 

@<partition@>=@#
static double  GenPartitions(const unsigned int m,  const unsigned int myInt,
                   const unsigned int PartitionSize,
                   unsigned int MinVal,
			  unsigned int MaxVal,
			  std::vector< std::pair< double, std::vector<unsigned int> > > & v_l_k,
			  std::vector<double>& lrates_sorting)
{

@#

  double lrate {} ;
  double sumrates {} ;
  
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

@#

        assert( static_cast<unsigned int>(std::accumulate( partition.begin(), partition.end(), 0)) == myInt);


   /* \newline \aths{\S~\ref{sec:lambdanks}} */

        lrate =  lambdanks( static_cast<double>(m), partition) ;

	assert( lrate >= 0) ;
	
	v_l_k.push_back( std::make_pair( lrate, partition) ) ;

	lrates_sorting.push_back( lrate );

	sumrates += lrate ; 
	

    } while (idx_Dec <= idx_Last);

    assert( sumrates >= 0) ;


/* \newline \aths{return the sum of the rates \eqref{eq:1} for the
generated partitions} */
    return sumrates ;
}


@*1 {\bf generate all partitions  summing to a given number}.
\label{sec:allmergersumm}


generate all partitions summing to a given number |m| when |n| is the
number of blocks 


@<fixed-sum partitions@>=@# 
static double allmergers_sum_m( const unsigned int n, const unsigned int m, std::vector< std::pair< double, std::vector<unsigned int> > >& v__l_k,  std::vector<double>&  v_lrates_sort )
{


@#

  /* \newline \aths{ |n| is number of blocks ;
      all partitions summing to $m \le  n$} 
   */
  const std::vector<unsigned int> v__m {m};
  /* \newline \aths{\S~\ref{sec:lambdanks}} */
  double sumr = lambdanks( static_cast<double>(n), v__m ) ;

@#

  v__l_k.push_back( std::make_pair( sumr, v__m ) ) ;

  v_lrates_sort.push_back( sumr );

@#
  if( m > 3){
  /* \newline \aths{|s| is the size of the partitions summing to |m|} */ 
    for( unsigned int s = 2; s <= m/2; ++s){
  @#    
      assert(m > 2*(s-1) );
      /* \newline \aths{\S~\ref{sec:partsumm}} */
      sumr += GenPartitions(n,  m, s, 2, m - (2*(s-1)), v__l_k, v_lrates_sort ); }
  }

@#

  assert( sumr >= 0) ;
  return sumr ;
}




@*1 {\bf write partitions to file}. 
\label{sec:ratesmergersfile}

write partitions  for given number of blocks |n|  to file


@<file partitions@>=@#
static void ratesmergersfile( const unsigned int n, const std::vector<unsigned int>& v__indx, const std::vector< std::pair< double, std::vector<unsigned int> > > & vlk, const double s,   std::vector< std::vector< double > > & a__cmf )
{

@#

  assert( s > 0);
  double cmf {} ;
  std::ofstream f ;
  f.open("gg_" + std::to_string(n) + "_.txt", std::ios::app);

@#

  a__cmf[n].clear() ;

@#

  for( const auto &i : v__indx)
    {

@#

      cmf += (vlk[ i ].first) / s ;
      assert( cmf >= 0) ;

@#

      a__cmf[n].push_back(cmf) ;
      assert((vlk[ i ].second).size() > 0 );

@#


      for( const auto &x : vlk[ i ].second){
	f << x << ' ' ;}
      f << '\n' ; }

@#

  f.close();

@#

  assert( abs(cmf - 1.) < 0.999999 );
}



@*1 {\bf generate all partitions when $n$ blocks}.
\label{sec:allmergerswhenblocks}



generate all partitions when given $n$ blocks 


@<generate all partitions when $n$ blocks@>=@#
static void allmergers_when_n_blocks( const unsigned int n, std::vector<double> & v__lambdan, std::vector< std::vector< double > > & a__cmf )
{

@#

  std::vector< std::pair< double, std::vector<unsigned int> > > vlk {} ;
  std::vector< double > ratetosort {} ;
  ratetosort.clear() ;

  double lambdan {} ;
  vlk.clear() ;
  assert( n > 1);

@#

  for( unsigned int k = 2 ; k  <= n ; ++k){
     /* \newline \aths{ \S~\ref{sec:allmergersumm};   the partition sums to |k|; the number of blocks is |n|} */
    lambdan += allmergers_sum_m(n, k, vlk, ratetosort); }
  /* \newline \aths{ record the total rate when |n| blocks} */
  /* \newline \aths{ use for sampling time} */
  assert( lambdan > 0);
  v__lambdan[n] = lambdan ;

@#
  
  std::vector<unsigned int> indx (ratetosort.size());

@#

  std::iota( indx.begin(), indx.end(), 0);


@#

  std::stable_sort( indx.begin(), indx.end(), [&ratetosort](const unsigned int x, const unsigned int y){return ratetosort[x] > ratetosort[y];});

 
 /* \newline \aths{\S~\ref{sec:ratesmergersfile} } */
  ratesmergersfile(n, indx, vlk,  v__lambdan[n], a__cmf);
}



@*1 {\bf all partitions}. 
\label{sec:allpartitionsrates}




 generate all partitions and rates  \eqref{eq:1}  up to sample size


@<generate all partitions and rates@>=@#
static void all_partitions_rates( std::vector<double>& vlmn, std::vector< std::vector<double> > & acmf  )
{
   
  for( unsigned int tmpn = 2; tmpn <= SAMPLE_SIZE; ++tmpn )
    {
      /* \newline \aths{\S~\ref{sec:allmergerswhenblocks}} */
      allmergers_when_n_blocks( tmpn, vlmn, acmf );
    }
}


@*1 {\bf sample a partition}.
\label{sec:samplemerger}

sample a partition given  number of blocks


@<sample partition@>=@#
static unsigned int samplemerger(const unsigned int n, const std::vector< double > & v__cmf )
{

@#

/* \newline \aths{|n| the number of blocks; |v__cmf| the cumulative
mass function based on \eqref{eq:1} } */


  unsigned int j {} ;
  const double u = gsl_rng_uniform( rngtype);

@#
  
  while( u > v__cmf[j]){ ++j; } 

@#
  
  return j ;
}



@*1 {\bf read partition from file}.
\label{sec:readmergersizes}


read a partition from file given the line indicator sampled using \S~\ref{sec:samplemerger}


@<fetch partition from file@>=@# 
static void readmergersizes(const unsigned int n,  const unsigned int j, std::vector<unsigned int> & v__mergers )
{

@#

  std::ifstream f("gg_" + std::to_string(n) + "_.txt") ;
  std::string line {} ;

@#

  v__mergers.clear(); 


@#

  for( unsigned int i = 0; std::getline(f, line) && i < j; ++i ){
    if(i >= j-1 ){

@#

      std::stringstream ss(line) ;

@#

      v__mergers = std::vector<unsigned int>( std::istream_iterator<unsigned int>(ss), {} );
    } }

@#

  assert(v__mergers.size() > 0) ;

@#

  assert( std::all_of(v__mergers.cbegin(), v__mergers.cend(), [](const auto &x){ return x > 1; }) );

@#

  f.close();
}



@*1 {\bf sample a box}.
\label{sec:samplebox}


sample one of four  boxes uniformly; blocks in the same box will be merged


@<pick a box@>=@#
unsigned int sample_box()
{

@#

  const double u = gsl_rng_uniform(rngtype) ;

@#

  return (u < 0.25 ? 0 : ( u < .5 ? 1 : ( u < .75 ? 2 : 3 ) ) ) ;
}



@*1 {\bf split one partition element}.
\label{sec:splitblocks}


split the blocks in one element   of  a partition using
\S~\ref{sec:samplebox}; we sample  $k$ iid random variables $I_{1},\ldots,
I_{k}$  where  $\prb{I_{i} = j} = 1/4$ for $j = 0,1,2,3$   


@<split partition element into boxes@>=@#
void split_blocks(const unsigned int k, std::vector<unsigned int>& split_partition )
{

@#
 
  /* \newline \aths{|k| is the given  size of the given partition  element} */ 

  std::vector<unsigned int> boxes (4);

@#

 /* \newline \aths{split the |k| blocks into boxes; $|boxes[j]|  = \sum_{i=1}^{|k| }\one{I_{i} = j}   $ } */
  for( unsigned int j = 0; j < k; ++j)
    {
    /* \newline \aths{\S~\ref{sec:samplebox}} */
      ++boxes[ sample_box() ];
    }

@#
    assert( static_cast<unsigned int>( std::accumulate(boxes.cbegin(), boxes.cend(), 0))  == k);

  /* \newline  \aths{ count how many boxes have at least 2 blocks} */
  const auto t = static_cast<std::size_t>( std::count_if( boxes.cbegin(), boxes.cend(), [](const auto b){ return b > 1;}));

@#

  if(t > 0){
  /* \newline \aths{append  box sizes with at least 2 blocks to |split_partition|} */
  std::copy_if(boxes.cbegin(), boxes.cend(), std::back_inserter(split_partition), [](const auto b){return b > 1;});
  } 
}



@*1 {\bf split a given partition into boxes}.
\label{sec:splitgroups}


split a partition into boxes; given a partition $(k_{1}, \ldots,
k_{r})$ split each of $k_{i}$ blocks into boxes and record the box
sizes with at least 2 blocks 


@<box a given partition@>=@#
void split_groups( const std::vector< unsigned int>& partition, std::vector<unsigned int>& split_partition )
{

@#

  split_partition.clear();

@#

  for( const auto &k:partition){
  /* \newline \aths{\S~\ref{sec:splitblocks}} */
    split_blocks( k,  split_partition ) ;
    }


@#

  assert( std::accumulate(split_partition.cbegin(), split_partition.cend(),0) <= std::accumulate(partition.cbegin(), partition.cend(),0));
}



@*1 {\bf update lengths $\ell_{i}$}.
\label{sec:updatelengths}


update realisation $\ell_{i}$ of $L_{i}(n)$; |t| the interval time
until next merger during which have tree configuration (block sizes)
|tree|  

@<given interval time update lengths@>=@#
static void update_lengths( const std::vector<unsigned int>& tree, std::vector<double>& v_l, const double t )
{

 @#

  for( const auto &b : tree){
    v_l[0] += t;
    v_l[b] += t; }
}



@*1 {\bf update $r_{i}$}.
\label{sec:updateri}

update approximations $r_{i}$ of $\EE{R_{i}(n)}$ given realised
lengths |v_l|   $\ell_{1},\ldots, \ell_{n-1}$



@<add to $r_{i}$@>=@#
static void update_ri( std::vector<double>& v_ri, const  std::vector<double>& v_l )
{

@#

  assert( v_l[0] > 0) ;
  const double d = v_l[0];


@#

   std::transform( v_l.begin(), v_l.end(), v_ri.begin(), v_ri.begin(), [&d](const auto &x, const auto &y){return y + (x/d);});
}




@*1 {\bf update tree}.
\label{sec:updatetree}


given merger size(s)   update  tree by merging blocks 


@<merge blocks and update tree@>=@#
static void update_tree( std::vector<unsigned int>& tree, const std::vector<unsigned int>& mergers )
{

@#

  assert(  mergers.size() > 0 );

@#

  assert( static_cast<std::size_t>(std::accumulate(mergers.cbegin(), mergers.cend(), 0)) <= tree.size() ) ;


@#

  std::shuffle( tree.begin(), tree.end(), rng) ;

@#

  /* \newline \aths{ |newblocks| record the size of the new blocks (obtained
  by merging blocks according to |mergers|  in the current tree)}  */
  std::vector<unsigned int> newblocks (mergers.size());

@#

  std::size_t j {} ;

@#

  for( const auto &m: mergers ){
  @#
    newblocks[j] = std::accumulate( std::crbegin(tree), std::crbegin(tree) + m, 0);

@#

    tree.resize( tree.size() - m ) ;

@#

    ++j ;
  }
  tree.reserve( tree.size() + newblocks.size() );

@#

  tree.insert( tree.end(), newblocks.cbegin(), newblocks.cend()  );
}



@*1 {\bf until a merger occurs}.
\label{sec:untilmerger}


the goal here is to 
\begin{enumerate}
\item sample time and  partitions and box sizes until have at least one box with at
least 2 blocks;
\item titions and box sizes until have at least one box with at
least 2 blocks; then shuffle the tree (the current block sizes)
\item  merge
blocks by summing block sizes  of the rightmost blocks of the shuffled
tree given merger sizes
\item remove blocks that have merged from the tree
\item append new blocks from the merged blocks  to the tree 
\end{enumerate}




@<sample time and boxes until merger@>=@#
static double until_merger(const std::size_t current_number_blocks,   const std::vector<double>& v_lambdan, const std::vector<double>& v__cmf, std::vector<unsigned int>& v_merger_sizes)
{

@#

  std::vector<unsigned int> v_partition  (SAMPLE_SIZE/2);


@#

  v_partition.reserve(SAMPLE_SIZE/2) ;


@#

  unsigned int lina {} ;
  double t {};

  v_merger_sizes.clear() ;

@#

  while( std::all_of(v_merger_sizes.cbegin(),  v_merger_sizes.cend(), [](const auto m){ return m < 2;}))
    {
      t += gsl_ran_exponential(rngtype, 1./v_lambdan[current_number_blocks] );
      /* \newline \aths{\S~\ref{sec:samplemerger}} */
      lina = samplemerger(current_number_blocks,  v__cmf);
       /* \newline \aths{\S~\ref{sec:readmergersizes}} */
      readmergersizes( current_number_blocks, 1 + lina, v_partition) ;
        /* \newline \aths{\S~\ref{sec:splitgroups}} */
       split_groups( v_partition, v_merger_sizes) ;
    }

@#

  assert( v_merger_sizes.size() > 0) ;

@#

  assert( std::all_of(v_merger_sizes.cbegin(),  v_merger_sizes.cend(), [](const auto m){ return m > 1;}));


@#
  return t ;
}



@*1 {\bf one experiment}.
\label{sec:oneexperiment}


generate one realisation of  $\ell_{1},\ldots, \ell_{n-1}$ starting
with tree configuration $(1,\ldots, 1)$ (all blocks of size 1)  and
ending with  one block of size $n$

@<generate one realisation of $\ell_{1},\ldots, \ell_{n-1}$@>=@#
static void one_experiment( std::vector<double>& v_l, std::vector<double>& v_r, const std::vector<double>& v__lambdan, const std::vector< std::vector< double > >&  a__cmf)
{

@#

  std::vector<unsigned int> tree (SAMPLE_SIZE, 1) ;


@#

  tree.reserve(SAMPLE_SIZE) ;


@#

  std::fill(v_l.begin(), v_l.end(),0);

@#
  
  double t {};

@#

  std::vector<unsigned int> v__merger_sizes {};

@#

  v__merger_sizes.clear() ;
  v__merger_sizes.reserve(2*SAMPLE_SIZE) ;


@#

  while( tree.size() > 1)
    {
    /* \newline \aths{\S~\ref{sec:untilmerger}} */
      t = until_merger(tree.size(), v__lambdan, a__cmf[tree.size()], v__merger_sizes );  
     /* \newline \aths{\S~\ref{sec:updatelengths}} */
      update_lengths(tree, v_l, t);
      /* \newline \aths{\S~\ref{sec:updatetree}} */
      update_tree(tree, v__merger_sizes);
    }
  assert( tree.back() == SAMPLE_SIZE);

/* \newline \aths{\S~\ref{sec:updateri}} */
  update_ri( v_r, v_l);
}



@*1 {\bf approximate  $\EE{R_{i}(n)}$ }. 
\label{sec:approximateeri}


 approximate $\EE{R_{i}(n)}$ for a given number of |EXPERIMENTS| 

@<go ahead -- approximate $\EE{R_{i}(n)}$@>=@#
static void approximate_eri()
{

@#

  std::vector<double> vri (SAMPLE_SIZE) ;

@#

vri.reserve(SAMPLE_SIZE) ;

@#
  std::vector<double> v__l (SAMPLE_SIZE);

@#

v__l.reserve(SAMPLE_SIZE) ;


@#

  std::vector<double> v__lambdan (SAMPLE_SIZE + 1) ;


@#

v__lambdan.reserve(SAMPLE_SIZE + 1) ;



@#

  std::vector< std::vector< double > > a__cmfs (SAMPLE_SIZE + 1, std::vector<double> {} ) ;


@#

 /* \newline \aths{\S~\ref{sec:allpartitionsrates}} */
  all_partitions_rates(v__lambdan, a__cmfs );

@#

  int r = EXPERIMENTS + 1 ;

@#
  while( --r > 0)
    {
    /* \newline \aths{\S~\ref{sec:oneexperiment}} */
      one_experiment(v__l,  vri, v__lambdan,  a__cmfs);
    }

@#

  for(const auto &z: vri){std::cout << z << '\n';}
}




@*1 {\bf main}.
\label{sec:main}

the |main|  module


@C


/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:descfact}} */
@<compute descending factorial@>@#
/* \newline \aths{\S~\ref{sec:power}} */
@<checked power function@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<compute  $\lambda_{n;k_{1},\ldots, k_{r};s}$ \eqref{eq:1}@>@#
/* \newline \aths{\S~\ref{sec:partsumm}}*/
@<partition@>@#
/* \newline \aths{\S~\ref{sec:allmergersumm}} */
@<fixed-sum partitions@>@# 
/* \newline \aths{\S~\ref{sec:ratesmergersfile}} */
@<file partitions@>@#
/* \newline \aths{\S~\ref{sec:allmergerswhenblocks}} */
@<generate all partitions when $n$ blocks@>@#
/* \newline \aths{\S~\ref{sec:allpartitionsrates}} */
@<generate all partitions and rates@>@#
/* \newline \aths{\S~\ref{sec:samplemerger}} */
@<sample partition@>@#
/* \newline \aths{\S~\ref{sec:readmergersizes}} */
@<fetch partition from file@>@# 
/* \newline \aths{\S~\ref{sec:samplebox}} */
@<pick a box@>@#
/* \newline \aths{\S~\ref{sec:splitblocks}} */
@<split partition element into boxes@>@#
/* \newline \aths{\S~\ref{sec:splitgroups}} */
@<box a given partition@>@#
/* \newline \aths{\S~\ref{sec:updatelengths}} */
@<given interval time update lengths@>@#
/* \newline \aths{\S~\ref{sec:updateri}} */
@<add to $r_{i}$@>@#
/* \newline \aths{\S~\ref{sec:updatetree}} */
@<merge blocks and update tree@>@#
/* \newline \aths{\S~\ref{sec:untilmerger}} */
@<sample time and boxes until merger@>@#
/* \newline \aths{\S~\ref{sec:oneexperiment}} */
@<generate one realisation of $\ell_{1},\ldots, \ell_{n-1}$@>@#
/* \newline \aths{\S~\ref{sec:approximateeri}} */
@<go ahead -- approximate $\EE{R_{i}(n)}$@>@#


int main(int argc, const  char * argv[])
{


/* \newline \aths{\S~\ref{sec:rngs}} */
 setup_rng( static_cast<std::size_t>(atoi(argv[1])));


 /* \newline \aths{\S~\ref{sec:approximateeri}} */
   approximate_eri() ;

  gsl_rng_free(rngtype);

  return 0 ;
}




@* {\bf conclusions}.
\label{sec:concl}


We approximate  $\EE{R_{i}(n)}$ for the
$\Omega$-$\delta_{0}$-Poisson-Dirichlet$(\alpha,0)$ coalescent with
$0< \alpha < 1$. Figure~\ref{fig:graph} records an example. 




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
        \addplot table[col sep=comma] {forplottingfile2};
        \addlegendentry{$\alpha = 0.01$, $c = 10^2$};
        \addlegendentry{$\alpha = 0.99$, $c = 10^0$};
       \end{axis}
       \end{tikzpicture}
       \caption{\it  An example approximation of $\EE{R_{i}(n)}$
graphed as logits  as a function of $\log(i/n) - \log(1 - i/n)$ for 
$i = 1,2,\ldots, n-1$ where $n$ is sample size }
      \label{fig:graph}
       \end{SCfigure}







@* {\bf bibliography}.
\label{sec:bib}


\bibliographystyle{plainnat}


\begin{thebibliography}{99}

\bibitem[D2024]{D2024}
 Diamantidis,  Dimitrios and Fan,  Wai-Tong (Louis) and Birkner,
 Matthias and Wakeley,  John.  {Bursts of coalescence within
 population pedigrees whenever big families occur}.  Genetics Volume 227,  February
  2024.
 \\
 \url{https://dx.doi.org/10.1093/genetics/iyae030}.


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
