%$ NAFN=timechanged_xindbetacoal_using_written_mergers
%$ ctangle $NAFN
%$ g++ -std=c++26 -m64 -march=native -O3 -x c++ timechanged_xindbetacoal_using_written_mergers.c -lm -lgsl -lgslcblas
%$ rm -f tmp_runs
%$ for i in $(seq 4); do echo "./a.out "  $(shuf -i 39393-282929191 -n1) "> resout"$i >> tmp_runs; done
%$ parallel --gnu -j4 :::: ./tmp_runs
%$ paste resout* -d, | sed '1d' | awk -F, '{M=1e4;s=0; for (i=1;i<=NF;i++) {s+=$i} print log(s/M) - log(1-(s/M))}' > logitresout
%$ seq 199 | awk '{n=200;print log($1/n) - log(1 - ($1/n))}' > nlogits
%$ paste -d',' nlogits logitresout > forplottingfile1
%$ cweave $NAFN
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
%$ lualatex $NAFN
%$ lualatex $NAFN
%% use write_upto_4fold_merger_sizes.cpp to generate the merger size
%% files Q_{m}_.txt
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
\usepackage{bbm}
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
\title{Gene genealogies in diploid populations evolving according to
sweepstakes reproduction \\ --- approximating $\EE{R_i(n)}$ for the
 time-changed $\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent}
 \author{Bjarki Eldon
 \footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{ compiled @@
{\DTMcurrenttime} on  {\today}  \\
\input{innleggctangle} \\
 \input{innleggcpp} \\ kernel  \input{innleggop} \\  \input{innleggbash} \\
GSL \input{innlegggsl} \\ \input{innleggcweave} \\  {\LaTeX}
\input{innlegglualatex}  \\  written using  \input{innleggemacs}}
\orcidlink{0000-0001-9354-2391}}


\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

\begin{abstract}
Let $\# A$ denote the number of elements in a finite set $A$.  For a
given coalescent $\set{\xi^{n}}$ write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :
\#\xi = i }dt$ and $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $
where $\tau(n) \equiv \inf \set{t \ge 0 : \#\xi^{n}(t) = 1}$. Then
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Write
$R_{i}(n) \equiv L_{i}(n)/L(n)$ for $i = 1,2,\ldots, n-1$.  With this
C++ simulation code we approximate $\EE{R_{i}(n)}$ when
$\set{\xi^{n}(G)} \equiv \set{\xi^{n}(G(t)) : t \ge 0 } $ is the
time-changed continuous-time  $\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$
coalescent with $0 < \alpha < 2$ and time-change function $G(t)$.
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
\label{sec:comp}

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


Write $\IN\equiv \set{1,2,\ldots}$, $[n] = \set{1,2,\ldots,n}$,
$]n] \equiv \set{2,3,\ldots,n}$ for
all $n\in \IN$. 
The $\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent is
a continuous-time $\Xi$-coalescent with  transition rates
\begin{equation}
\label{eq:1}
\begin{split}
\lambda_{n;k_{1},\ldots,k_{r};s} & =   \one{r=1,k_{1}=2}\frac{C_{\kappa}}{C_{\alpha,\gamma} } + \\
& +  \frac{c \alpha }{C_{\alpha,\gamma} \mathbbm{m}^{\alpha}}\sum_{\ell=0}^{s\wedge (4-r)} \binom{s}{\ell}\frac{(4)_{r+\ell} }{4^{k+\ell}}\int_{0}^{1}\one{0 < t \le \gamma }t^{k+\ell - \alpha - 1} {(1-t)}^{n + \alpha -k-\ell - 1}dt
\end{split}
\end{equation}
where  $0 < \gamma \le 1$, $B(p,a,b) \equiv  \int_{0}^{1}\one{0 < t\le p} t^{a-1}{(1-t)}^{b-1}dt$, and
\begin{displaymath}
\begin{split}
C_{\kappa} & =  \frac{2}{4\mathbbm{m}^{2} }\svigi{  \one{\kappa = 2} +  \one{\kappa > 2} \frac{c_{\kappa}}{2^{\kappa}(\kappa - 2)(\kappa - 1)}} \\
C_{\kappa, \alpha, \gamma} & = C_{\kappa} +  \frac{c\alpha}{4\mathbbm{m} ^{\alpha} }B(\gamma,2-\alpha,\alpha) \\
\mathbbm{m} & =  \one{\kappa = 2}\svigi{2\frac{\pi^{2}}{3} - 3}   +  \one{\kappa > 2}\svigi{1 +  2^{\kappa - 1}\frac{\kappa}{\kappa - 1} } 
\end{split}
\end{displaymath}
In \eqref{eq:1}  $n\ge 2$,  $k_{1},\ldots,k_{r}\in ]n]$  and $2 \le
\sum_{i}k_{i} \le n$  for all $r\in
[4]$, $s \equiv  n - \sum_{i}k_{i}$.   

Here we consider a time-changed version of the
$\Omega$-$\delta_{0}$-Beta$(\gamma,2-\alpha,\alpha)$-coalescent
corresponding to exponential population growth where the time-change
function $G(t) = \int_{0}^{t}e^{\rho s}ds $ for some fixed $\rho \ge
0$.   


Let $\# A$ denote the number of elements in a finite set $A$.  For a
given coalescent $\set{\xi^{n}}$ write
$L_{i}(n) \equiv \int_{0}^{\tau(n)} \# \set{ \xi \in \xi^{n}(t) :
\#\xi = i }dt$ and $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $
where $\tau(n) \equiv \inf \set{t \ge 0 : \#\xi^{n}(t) = 1}$. Then
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$. Write
$R_{i}(n) \equiv L_{i}(n)/L(n)$ for $i = 1,2,\ldots, n-1$.  With this
C++ code we use 
simulations to approximate $\EE{R_{i}(n)}$


 The code follows in
\S~\ref{sec:includes}--\S~\ref{sec:main}, we conclude in \S~\ref{sec:concl}. Comments within the code in
\aths{this font and color}


@* {\bf code}.
\label{sec:code}


@*1 {\bf includes}.
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





@*1 {\bf constants}.
\label{sec:constants}


the global constants 


@<constants@>=@#
/*
\newline 
  \aths{Model 0: $1 \le \alpha < 2$ : epsilon = cN**($\alpha - 2$)(1{$(\kappa > 2)$} + {$(\kappa = 2)$}log N);}
  */
  /* \newline \aths{
Model 1 : $\alpha = 1$ : epsilon = constant \\
$0 < \alpha < 1$ : epsilon = constant * N * (cutoff)**(alpha - 1)}
*/
const double dalpha = 0.01;
const long double ALPHA = static_cast<long double>(dalpha) ;
const long double KAPPA = 2.0L ;
const long double RHO = 0.0L ;
const double dgamma = 0.1;
const long double C_C = 1.0L ;
/* \newline \aths{ when $|KAPPA| > 2$:}  
const long double MM = $1 +  2^{\kappa} * \kappa / (\kappa - 1)$;
*/
/* \newline \aths{when  $|KAPPA| = 2$:}  */ 
const long double MM = (2.0L * M_PIl*M_PIl/3.0L) - 3.0L ;
const long double CA = ((KAPPA + 2.0L) + (KAPPA*KAPPA))/2.0L;
const long double CB = powl(2.0L , KAPPA)*((KAPPA - 2.0L)*(KAPPA - 1.0L));
const long double CKAPPA = 2.0L*(KAPPA > 2.0L ? CA/CB : 1.0L)/(MM*MM);
const long double GAMMA = static_cast<long double>(dgamma) ;
const long double BETA = static_cast<long double>(gsl_sf_beta(2. - dalpha, dalpha)*(dgamma < 1. ? gsl_sf_beta_inc(2.-dalpha, dalpha, dgamma) : 1.));
const long double CKAG = (CKAPPA + ((ALPHA*C_C * powl(2.0L, ALPHA))*BETA/powl(MM,ALPHA)/4.0L));
/* \newline \aths{|SAMPLE_SIZE| is number of diploid individuals
sampled} */
const unsigned int SAMPLE_SIZE = 100 ;
const int EXPERIMENTS = 2500 ;



@*1 {\bf random number generators}.
\label{sec:rngs}


the random number generators


@<rngs@>=@#
 std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine} */
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


@*1 {\bf compute $e^{x}$ checking for under and overflow}.
\label{sec:veldi}

compute $e^{x}$ checking for under and overflow


\begin{Piton}
long double veldi( const long double v )
{
  feclearexcept(FE_ALL_EXCEPT);

  const long double svar =  expl( v ) ;

  return( fetestexcept( FE_OVERFLOW ) != 0 ? LDBL_MAX : ( fetestexcept(FE_UNDERFLOW) != 0 ? 0.0L : svar) ) ;
}
\end{Piton}

@<$e^{x}$ with checks@>=@#
long double veldi( const long double v )
{
  feclearexcept(FE_ALL_EXCEPT);

  const long double svar =  expl( v ) ;

  return( fetestexcept( FE_OVERFLOW ) != 0 ? LDBL_MAX : ( fetestexcept(FE_UNDERFLOW) != 0 ? 0.0L : svar) ) ;
}


@*1 {\bf number of collisions}.
\label{sec:numbercollisions}


get the factorial number of collisions
$\prod_{j=2}^{n}\svigi{\sum_{i}\one{k_{i} = j}}!$

@<numbercollisions@>=@#
static long double numbercollisions( const std::vector<unsigned int>& __k )
{

/* \newline \aths{  |__k| is in ascending order} */
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

  assert(l > DBL_EPSILON) ;

  return static_cast<long double>(l) ;
}



@*1 {\bf multinomial constant}.
\label{sec:multinomialconstant}



compute  the log of  the multinomial constant
\begin{displaymath}
\binom{m}{k_{1}\ldots k_{r}s} \frac{1}{ \prod_{j=2}^{n}\svigi{\sum_{i}\one{k_{i} = j}}!  }
\end{displaymath}



@<multinomial constant@>=@#
long double multinomialconstant( const unsigned int m, const std::vector<unsigned int>& v__k)
{
@#
  const long double s = lgammal( static_cast<long double>(1 +  m - std::accumulate( v__k.cbegin(), v__k.cend(), 0)) ) ;
  
@#

  const long double d =   static_cast<long double>(std::accumulate( v__k.cbegin(), v__k.cend(), 0.0L, [](long double a, const auto x){return a + lgammal(static_cast<long double>(x+1)); })) ;

/* \newline \aths{|numbercollisions| \S~\ref{sec:numbercollisions}} */
  return ( lgammal(static_cast<long double>(m + 1)) - d - s - logl(numbercollisions(v__k)) ) ;
}



@*1 {\bf the (incomplete) beta function}.
\label{sec:betafunc}



 return the logarithm of the (incomplete) beta function; when complete
 it  is  of course   $\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b)$ 

@<the beta function@>=@#
static long double betafunc (const double x, const double a, const double b)
{
  /* \newline \aths{ the GSL incomplete beta function is normalised by the complete beta function} */
  /* \newline \aths{ $0 < x <= 1$  is the cutoff point} */
  assert( x > 0.);
  assert( a > 0.);
  assert( b > 0.);
  /* \newline the standard way would be \aths{ | gsl_sf_beta_inc( a, b, x) * gsl_sf_beta( a, b) |} */

  const long  double f = static_cast<long double>( (x < 1. ? gsl_sf_hyperg_2F1(a + b, 1, a+1.,x) : 1.) );
  
  assert( LDBL_EPSILON < f);
  /* \newline  return $\one{x < 1}( \log f  +a\log x + b\log(1-x) - \log a) + \one{x = 1} (\log\Gamma(a) + \log\Gamma(b) - \log\Gamma(a+b))   $ */
  return( x < 1. ?  (logl( f ) +  ( static_cast<long double>( (a*log(x)) + (b * log(1-x)) - log(a) ) ) ) : (lgammal( static_cast<long double>(a)) + lgammal( static_cast<long double>(b) ) - lgammal( static_cast<long double>(a + b) ) ) ); 
}



@*1 {\bf read merger sizes}.
\label{sec:readmergersizes}


read in merger sizes summing to a given number; the merger sizes are
available in the files \\  {\tt Q\_<m>\_.txt}  \\ where {\sf m} is the
given number


@<readmergersizes@>=@#
static void readmergersizes(const unsigned int n,  const unsigned int j, std::vector<unsigned int> & v__mergers )
{
@#
  std::ifstream f("Q_" + std::to_string(n) + "_.txt") ;

@#

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


@*1 {\bf $\lambda_{n;k_{1},\ldots,k_{r};s}$}.
\label{sec:lambdanks}


compute  $\lambda_{n;k_{1},\ldots,k_{r};s}$ \eqref{eq:1}


@<$\lambda_{n;k_{1},\ldots,k_{r};s}$@>=@#
long double lambdanks( const unsigned int m, const std::vector<unsigned int>& v__k )
{

@#

/* \newline \aths{$r$ is number of simultaneous mergers} */
  const unsigned int  r = v__k.size();
  /* \newline \aths{$k = k_{1} + \cdots + k_{r}$ } */
  const unsigned int k = std::accumulate( v__k.cbegin(),  v__k.cend(),0);
  assert( k > 1) ;
  assert( k <= m);

@#

  const unsigned int s = m - k ;

@#

/* \newline \aths{$\log\binom{x}{y} = \log\Gamma(x+1) -
\log\Gamma(y+1) - \log\Gamma(x-y+1)$} */
  auto logchoose = [](const unsigned int x, const unsigned int y){ return (lgammal( static_cast<long double>(x+1)) - lgammal(static_cast<long double>(y + 1)) - lgammal(static_cast<long double>(x-y + 1)));};

@#

  long double l = 0.0L;

@#

/* \newline \aths{$(4)_{p} \equiv 4(4-1)\cdots (4-p+1) $ and
$(4)_{0}\equiv 1$} */
auto ff = [](const unsigned int p){
return (p > 0 ? static_cast<long double>(boost::math::falling_factorial(4.0, p)) : 1.0L);};

@#

  for( std::size_t ell = 0; ell <=  (s < 4-r ? s : 4-r); ++ell ){
    l += veldi(logchoose(s,ell)  +  ff(static_cast<long double>(r+ell))  + betafunc(GAMMA, static_cast<long double>(k + ell) - ALPHA, static_cast<long double>(m-k-ell)+ALPHA) - (logl(4.0L)*static_cast<long double>(k+ell)) ); }

@#

  l *= (ALPHA * C_C) * powl(2.0L, ALPHA) ;
  l /= (CKAG * powl( MM, ALPHA) );

  l += (k < 3 ?  CKAPPA / CKAG : 0.0L);
  assert( l > 0.0L ) ;
  
   return ( veldi( logl(l) +  multinomialconstant( m, v__k ) ) ) ;
}



@*1 {\bf rates for a given number of blocks}.
\label{sec:allrates_when_n_blocks}


read in all the possible merger sizes  when a given number of blocks and
compute the rate \eqref{eq:1}

@<rates when given number of blocks@>=@#
void allrates_when_n_blocks(const unsigned int n, std::vector<long double>& v__lambdan, std::vector< long double>& v__rates )
{
@#
  std::vector< unsigned int> v__k {} ;
  v__k.clear();
  std::string line {};
  std::ifstream f ("Q_" + std::to_string(n) + "_.txt");
  assert( f.is_open() );
  std::stringstream ss {};
  long double l {};
  v__rates.clear();
  @#
  while( std::getline(f, line) ){
    @#
    ss = std::stringstream(line) ;
    @#
    v__k = std::vector<unsigned int>( std::istream_iterator<unsigned int>(ss), {} );
    @#
    assert(v__k.size() > 0) ;
    @#
    assert( std::all_of(v__k.cbegin(), v__k.cend(), [](const auto &x){return x > 1; }) );
    /* \newline \aths{compute $\lambda_{n;k_{1},\ldots,k_{r};s}$ \eqref{eq:1} \S~\ref{sec:lambdanks}} */
    l = lambdanks(n, v__k);
    v__lambdan[n] += l ;
    v__rates.push_back(l); 
  }
  @#
  f.close();
}



@*1 {\bf order the rates in descending order}. 
\label{sec:order_rates}



order the rates in descending order for  generating the cmf for
sampling merger sizes


\begin{Piton}
static void order_rates(  const std::vector< long double > & v__rates_for_sorting, std::vector< unsigned int > & v__indx )
{
  assert( v__rates_for_sorting.size() > 0) ;
  v__indx.clear();
  v__indx.resize( v__rates_for_sorting.size() ) ;
  std::iota( v__indx.begin(), v__indx.end(), 0 ) ;
  std::stable_sort(v__indx.begin(), v__indx.end(), [&v__rates_for_sorting](const unsigned int x, const unsigned int y){return v__rates_for_sorting[x] > v__rates_for_sorting[y];});
}
\end{Piton}

@<order rates@>=@#
static void order_rates(  const std::vector< long double > & v__rates_for_sorting, std::vector< unsigned int > & v__indx )
{
@#
  assert( v__rates_for_sorting.size() > 0) ;
  v__indx.clear();
  v__indx.resize( v__rates_for_sorting.size() ) ;
  std::iota( v__indx.begin(), v__indx.end(), 0 ) ;
  std::stable_sort(v__indx.begin(), v__indx.end(), [&v__rates_for_sorting](const unsigned int x, const unsigned int y){return v__rates_for_sorting[x] > v__rates_for_sorting[y];});
}



@*1 {\bf generate cmf}.
\label{sec:generate_cmf}


 generate cmf for a given set of merger sizes

\begin{Piton}
static void generate_cmf_n_blocks( const long double lambdan,   const std::vector<long double>& v__rates, const std::vector<unsigned int> & v__indx, std::vector<long double>& v__cmf )
{
  v__cmf.clear();
  v__cmf.resize( v__rates.size() ); 
  v__cmf[0] = v__rates[v__indx[0]] / lambdan ;
  for( unsigned int i = 1; i < v__indx.size(); ++i)
    {
      v__cmf[i] = v__cmf[ i-1] + ( v__rates[ v__indx[i] ] / lambdan ); 
    }
  assert( fabsl(v__cmf.back() - 1.0L) < 1.0e-9L ) ;
}
\end{Piton}


@<generate cmf@>=@#
static void generate_cmf_n_blocks( const long double lambdan,   const std::vector<long double>& v__rates, const std::vector<unsigned int> & v__indx, std::vector<long double>& v__cmf )
{
/* \newline \aths{|lambdan| is the total merger rate} */
  v__cmf.clear();
  v__cmf.resize( v__rates.size() ); 
  v__cmf[0] = v__rates[v__indx[0]] / lambdan ;
  @#
  for( unsigned int i = 1; i < v__indx.size(); ++i)
    {
      v__cmf[i] = v__cmf[ i-1] + ( v__rates[ v__indx[i] ] / lambdan ); 
    }
  assert( fabsl(v__cmf.back() - 1.0L) < 1.0e-9L ) ;
}



@*1 {\bf compute all merger rates}. 
\label{sec:allrates}


compute the rate \eqref{eq:1}  for all possible mergers

@<all rates@>=@#
void allrates( std::vector<long double>& v__lambdan, std::vector< std::vector< unsigned int> > & a__indx, std::vector< std::vector< long double> > & a__cmfs)
{

  std::vector< std::vector< long double > > a__rates (SAMPLE_SIZE + SAMPLE_SIZE + 1, std::vector< long double> {}) ;

  /* \newline  \aths{|SAMPLE_SIZE| \S~\ref{sec:constants}  is number of diploid individuals sampled}  */
  for( unsigned int m = 2; m <= SAMPLE_SIZE + SAMPLE_SIZE; ++m)
    {
    /* \newline \aths{\S~\ref{sec:allrates_when_n_blocks}} */
      allrates_when_n_blocks(m, v__lambdan, a__rates[m]);
      assert( v__lambdan[m] > LDBL_EPSILON) ;
      /* \newline \aths{\S~\ref{sec:order_rates}} */
      order_rates(a__rates[m], a__indx[m]);
      /* \newline \aths{\S~\ref{sec:generate_cmf}} */
      generate_cmf_n_blocks(v__lambdan[m], a__rates[m], a__indx[m], a__cmfs[m]);
      a__rates[m].resize(0);
      @#
      std::vector<long double>().swap(a__rates[m]);
    }

/* \newline \aths{clean up}*/
  for( auto &v: a__rates){
    v.resize(0);
    std::vector<long double>().swap(v);}
  std::vector< std::vector< long double> >().swap(a__rates);
}



@*1 {\bf get time until next merger}.
\label{sec:newtime}


get the time
\begin{displaymath}
t =  \one{\rho > 0} \frac{1}\rho  \log\svigi{1 -\rho \log(1-U)  e^{-\rho s} / \lambda_{n}} + \one{\rho = 0} \text{Exp}(\lambda_{n}) 
\end{displaymath}
where $U$ is a standard random uniform and $s$ is the sum of the
previous times

\begin{Piton}
long double newtime(const long double lambdab, const long double oldtime)
{
  return (RHO > LDBL_EPSILON ? log1pl( -(RHO * expl(- RHO * oldtime) / lambdab) * logl(static_cast<long double>(gsl_rng_uniform_pos(rngtype))))/RHO : static_cast<long double>(gsl_ran_exponential(rngtype, 1./lambdab))) ;
}
\end{Piton}


@<new time@>=@#
long double newtime(const long double lambdab, const long double oldtime)
{
/* \newline \aths{|log1p(x)| is $\log(1+x)$} */
  return (RHO > LDBL_EPSILON ? log1pl( -(RHO * expl(- RHO * oldtime) / lambdab) * logl(static_cast<long double>(gsl_rng_uniform_pos(rngtype))))/RHO : static_cast<long double>(gsl_ran_exponential(rngtype, 1./lambdab))) ;
}



@*1 {\bf update lengths}.
\label{sec:update_lengths}


update functionals $L_{i}(n)$   given  current block sizes


@<update lengths@>=@#
void update_lengths( const std::vector<unsigned int>& v__tree,   std::vector<long double>& v__lengths, const long double t )
{
@# 
  for( unsigned int i = 0; i < v__lengths.size(); ++i)
    {
      v__lengths[0] += t ;
      assert(v__tree[i] > 0 ) ;
      assert(v__tree[i] < SAMPLE_SIZE + SAMPLE_SIZE) ;
      v__lengths[v__tree[i]] += t ;
    }
}


@*1 {\bf update block sizes}.
\label{sec:update_tree}

merge blocks and record the continuing ones

\begin{Piton}
void update_tree( std::vector<unsigned int> & v__tree, const std::vector<unsigned int>& merger_sizes )
{
  assert( merger_sizes.size() > 0) ;
  assert( static_cast<unsigned int>(std::accumulate( merger_sizes.cbegin(), merger_sizes.cend(),0)) <= v__tree.size() );
  std::shuffle( v__tree.begin(), v__tree.end(), rng );
  std::vector<unsigned int > newblocks {};
  newblocks.clear();
  newblocks.reserve( merger_sizes.size() ) ;
  assert( newblocks.size() < 1) ;
  for( const auto &s : merger_sizes){
    newblocks.push_back( std::accumulate( v__tree.crbegin(), v__tree.crbegin() + s, 0));
    v__tree.resize( v__tree.size() - s) ; }

  v__tree.insert( v__tree.end(), newblocks.cbegin(), newblocks.cend() );  
}
\end{Piton}


@<update block sizes@>=@#
void update_tree( std::vector<unsigned int> & v__tree, const std::vector<unsigned int>& merger_sizes )
{
@#
  assert( merger_sizes.size() > 0) ;
  @#
  assert( static_cast<unsigned int>(std::accumulate( merger_sizes.cbegin(), merger_sizes.cend(),0)) <= v__tree.size() );

@#

std::shuffle( v__tree.begin(), v__tree.end(), rng );
  std::vector<unsigned int > newblocks {};
  newblocks.clear();
  newblocks.reserve( merger_sizes.size() ) ;
  assert( newblocks.size() < 1) ;

@#

  for( const auto &s : merger_sizes){
    newblocks.push_back( std::accumulate( v__tree.crbegin(), v__tree.crbegin() + s, 0));
    v__tree.resize( v__tree.size() - s) ; }

  v__tree.insert( v__tree.end(), newblocks.cbegin(), newblocks.cend() );  
}


@*1 {\bf one experiment}. 
\label{sec:one_experiment}


\begin{Piton}
void one_experiment( const std::vector<long double>& v__lambdan,  std::vector<long double>&  v__lengths, const  std::vector< std::vector< long double > > & a__cmfs,  const std::vector< std::vector< unsigned int> > & a__indx, std::vector<long double>& v__ri )
{
  std::vector< unsigned int> v__tree( 2*SAMPLE_SIZE, 1);
  std::fill( v__lengths.begin(), v__lengths.end(), 0.0L) ;
  unsigned int merger_lina {} ;
  std::vector< unsigned int > merger_sizes {} ;
  long double otimi {};
  long double it {};


  auto mlina = [](const std::vector<long double>& f)
  {
    unsigned int j = 0;
    const long double u = static_cast<long double>(gsl_rng_uniform(rngtype));
    assert( u <= 1.0L );

  while( u > f[j]){++j;}

  return j ;};

  
  merger_sizes.reserve(4);
  while( v__tree.size() > 1)
    {
      it = newtime( v__lambdan[v__tree.size()],  otimi);
      otimi += it ;
      update_lengths(v__tree, v__lengths, it);
      merger_lina = mlina( a__cmfs[v__tree.size()]  ) ;
      readmergersizes( v__tree.size(), a__indx[v__tree.size()][merger_lina]+1,  merger_sizes);
      update_tree( v__tree, merger_sizes ) ;
    }
  assert( v__tree.back() == (2*SAMPLE_SIZE));

  assert( v__lengths[0] > LDBL_EPSILON) ;
  const long double d = v__lengths[0];
  std::transform( v__lengths.cbegin(), v__lengths.cend(), v__ri.begin(), v__ri.begin(), [&d](const auto x, const auto y){return y + (x/d);});
}
\end{Piton}


@<one experiment@>=@#
void one_experiment( const std::vector<long double>& v__lambdan,  std::vector<long double>&  v__lengths, const  std::vector< std::vector< long double > > & a__cmfs,  const std::vector< std::vector< unsigned int> > & a__indx, std::vector<long double>& v__ri )
{
@#
  std::vector< unsigned int> v__tree( 2*SAMPLE_SIZE, 1);
  std::fill( v__lengths.begin(), v__lengths.end(), 0.0L) ;
  unsigned int merger_lina {} ;
  std::vector< unsigned int > merger_sizes {} ;
  long double otimi {};
  long double it {};


  auto mlina = [](const std::vector<long double>& f)
  {
    unsigned int j = 0;
    const long double u = static_cast<long double>(gsl_rng_uniform(rngtype));
    assert( u <= 1.0L );

  while( u > f[j]){++j;}

  return j ;};

  
  merger_sizes.reserve(4);
  @#
  while( v__tree.size() > 1)
    {
    /* \newline \aths{\S~\ref{sec:newtime}} */
      it = newtime( v__lambdan[v__tree.size()],  otimi);
      otimi += it ;
      /* \newline \aths{\S~\ref{sec:update_lengths}} */
      update_lengths(v__tree, v__lengths, it);
      merger_lina = mlina( a__cmfs[v__tree.size()]  ) ;
      /* \newline \aths{\S~\ref{sec:readmergersizes}} */
      readmergersizes( v__tree.size(), a__indx[v__tree.size()][merger_lina]+1,  merger_sizes);
      /* \newline \aths{\S~\ref{sec:update_tree}} */
       update_tree( v__tree, merger_sizes ) ;
    }
  assert( v__tree.back() == (2*SAMPLE_SIZE));
@#
  assert( v__lengths[0] > LDBL_EPSILON) ;

@#

  const long double d = v__lengths[0];

@#

  std::transform( v__lengths.cbegin(), v__lengths.cend(), v__ri.begin(), v__ri.begin(), [&d](const auto x, const auto y){return y + (x/d);});
}



@*1 {\bf approximate}. 
\label{sec:approximate}


\begin{Piton}
void approximate()
{
  std::vector<long double> v__lambdan ((2*SAMPLE_SIZE) + 1) ;

  std::vector< std::vector< long double > > a__cmfs ( (2*SAMPLE_SIZE) + 1, std::vector< long double > {}) ;
  std::vector< std::vector< unsigned int > > a__indx ( (2*SAMPLE_SIZE) + 1, std::vector< unsigned int > {}) ;
  std::vector<long double> v__lengths (2*SAMPLE_SIZE) ;
  std::vector< long double> v__ri (2*SAMPLE_SIZE) ;

  allrates(v__lambdan, a__indx, a__cmfs);
  
  int r = EXPERIMENTS + 1;

  while( --r > 0)
    {one_experiment(v__lambdan, v__lengths, a__cmfs,  a__indx, v__ri) ; }

  for( const auto &z:v__ri)
    {std::cout << z << '\n';}
}
\end{Piton}


@<go ahead -- approximate $\EE{R_{i}(n)}$@>=@#
void approximate()
{
@#
  std::vector<long double> v__lambdan ((2*SAMPLE_SIZE) + 1) ;

  std::vector< std::vector< long double > > a__cmfs ( (2*SAMPLE_SIZE) + 1, std::vector< long double > {}) ;
  std::vector< std::vector< unsigned int > > a__indx ( (2*SAMPLE_SIZE) + 1, std::vector< unsigned int > {}) ;
  std::vector<long double> v__lengths (2*SAMPLE_SIZE) ;
  std::vector< long double> v__ri (2*SAMPLE_SIZE) ;

/*\newline \aths{\S~\ref{sec:allrates}} */
  allrates(v__lambdan, a__indx, a__cmfs);
  
  int r = EXPERIMENTS + 1;
@#

  while( --r > 0)
    {/* \newline \aths{\S~\ref{sec:one_experiment}} */
    one_experiment(v__lambdan, v__lengths, a__cmfs,  a__indx, v__ri) ; }

  for( const auto &z:v__ri)
    {std::cout << z << '\n';}
}






@*1 {\bf main}.
\label{sec:main}

the |main| module

@C

/* \newline \aths{\S~\ref{sec:includes}} */
@<includes@>@#
/* \newline \aths{\S~\ref{sec:constants}} */
@<constants@>@#
/* \newline \aths{\S~\ref{sec:rngs}} */
@<rngs@>@#
/* \newline \aths{\S~\ref{sec:veldi}} */
@<$e^{x}$ with checks@>@#
/* \newline \aths{\S~\ref{sec:numbercollisions}} */
@<numbercollisions@>@#
/* \newline \aths{\S~\ref{sec:multinomialconstant}} */
@<multinomial constant@>@#
/* \newline \aths{\S~\ref{sec:betafunc}} */
@<the beta function@>@#
/* \newline \aths{\S~\ref{sec:readmergersizes}} */
@<readmergersizes@>@#
/* \newline \aths{\S~\ref{sec:lambdanks}} */
@<$\lambda_{n;k_{1},\ldots,k_{r};s}$@>@#
/* \newline \aths{\S~\ref{sec:allrates_when_n_blocks}} */
@<rates when given number of blocks@>@#
/* \newline \aths{\S~\ref{sec:order_rates}} */
@<order rates@>@#
/* \newline \aths{\S~\ref{sec:generate_cmf}} */
@<generate cmf@>@#
/* \newline \aths{\S~\ref{sec:allrates}} */
@<all rates@>@#
/* \newline \aths{\S~\ref{sec:newtime}} */
@<new time@>@#
/* \newline \aths{\S~\ref{sec:update_lengths}} */
@<update lengths@>@#
/* \newline \aths{\S~\ref{sec:update_tree}} */
@<update block sizes@>@#
/* \newline \aths{\S~\ref{sec:one_experiment}} */
@<one experiment@>@#
/* \newline \aths{\S~\ref{sec:approximate}} */
@<go ahead -- approximate $\EE{R_{i}(n)}$@>@#

int main(int argc, const char * argv[])
{


/* \newline \aths{\S~\ref{sec:rngs}} */
  setup_rng( static_cast< std::size_t>( atoi(argv[1])) ) ;

/* \newline \aths{\S~\ref{sec:approximate}} */
approximate();


/* \newline \aths{|rngtype| \S~\ref{sec:rngs}} */
gsl_rng_free (rngtype);
return GSL_SUCCESS;
}


@* {\bf conclusion and references}.
\label{sec:concl}


Figure~\ref{fig:graph} records an example approximation of
$\EE{R_{i}(n)}$ given the  parameter values in \S~\ref{sec:constants}





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
        %\addlegendentry{$\alpha = 0.05$}
        %\addlegendentry{$\alpha = 1$}
        %\addlegendentry{$\alpha = 1.5$}
       \end{axis}
       \end{tikzpicture}
       \caption{\textcolor{violet}{ \it  An example approximation of $\EE{R_{i}(n)}$
 for the  given  parameter values and graphed as logits against
 $\log(i/n) - \log(1 - i/n)$ for $i = 1,2,\ldots, n-1$ where $n$ is
 sample size}}
       \label{fig:graph}
       \end{SCfigure}








@* {\bf bibliography}.
\label{sec:bib}


\bibliographystyle{plain}


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
