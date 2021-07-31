(* ::Package:: *)

(* ::Text:: *)
(*This Mathematica notebook refers to the following paper:*)
(**)
(*I. Palaia, A. Paraschiv, V. Debets, C. Storm, A. \[CapitalSHacek]ari\[CAcute]*)
(*Durotaxis of passive nanoparticles on elastic membranes*)
(*bioR\[Chi]iv (2021)*)
(*https://doi.org/10.1101/2021.04.01.438065*)
(**)
(*It contains code to compute the free energy of a passive nanoparticle adhering to a fluctuating, bendable membrane. *)
(*First, adhesion energy, bending energy, global stretching, and fluctuation entropy are defined as a function of wrapped area (see Methods section of the paper). *)
(*Then, the constrained free energy resulting from the sum of these terms is minimised, and the equilibrium free energy (see Fig. 3) is output to a file.*)
(**)
(*Ivan Palaia, 12Jul2021, i.palaia@ucl.ac.uk*)
(**)


(* ::Section:: *)
(*F minimization*)


(* ::Text:: *)
(*Notation:*)
(**)
(*\[CapitalSigma] is wrapped area, in units of \[Sigma]^2 where \[Sigma] is the unit of length (e.g. in simulations).*)
(*\[Kappa] is bending rigidity, in units of thermal energy Subscript[k, B]T.*)
(*\[Tau] is surface tension, in units of Subscript[k, B]T/\[Sigma]^2.*)
(*a is the side of a square whose area equals the area of a small rhombus, defined by the membrane mesh in simulations. In the paper, a = n^(-1/2). It is in units of \[Sigma].*)
(*\[Epsilon] is the adhesion energy per bead, in units of Subscript[k, B]T. In other words, \[Epsilon] and a are such that the adhesion energy per surface is \[Epsilon]/a^2. See Eq. (6).*)
(*R is the radius of the adhered nanoparticle, in units of \[Sigma].*)
(**)
(*\[Delta]h2Ten is the squared amplitude of fluctuations divided by Subscript[k, B]T, as defined by Eq. (11). It is in units of \[Sigma]^2/(Subscript[k, B]T).*)
(*F is the constrained free energy, in units of Subscript[k, B]T.*)
(*Ebend is the bending energy, in units of Subscript[k, B]T.*)
(*Eadh is the adhesion energy, in units of Subscript[k, B]T.*)
(*Esurf is the energy associated with surface tension, in units of Subscript[k, B]T.*)
(*Sm is the entropic term of the free energy -TS, in units of Subscript[k, B]T.*)
(**)
(**)


(* ::Input::Initialization:: *)
\[Delta]h2Ten[\[CapitalSigma]_,a_,\[Epsilon]_,\[Kappa]_,\[Tau]_]:=((-ArcTan[(8 \[Pi]^2 \[Kappa]+\[CapitalSigma] \[Tau])/(\[CapitalSigma] Sqrt[(248 2^(2/3) \[Epsilon] \[Kappa])/a^4-\[Tau]^2])]+ArcTan[(8 \[Pi]^2 \[Kappa]+a^2 \[Tau])/Sqrt[248 2^(2/3) \[Epsilon] \[Kappa]-a^4 \[Tau]^2]])/(2 \[Pi] Sqrt[(248 2^(2/3) \[Epsilon] \[Kappa])/a^4-\[Tau]^2]))    


(* ::Input::Initialization:: *)
F[\[CapitalSigma]_,a_,\[Epsilon]_,\[Kappa]_,\[Tau]_,R_]:=(2\[Kappa])/R^2 \[CapitalSigma] - \[CapitalSigma] \[Epsilon]/a^2+\[Tau] \[CapitalSigma]^2/(4\[Pi] R^2)-\[CapitalSigma]/a^2  Log[Sqrt[\[Delta]h2Ten[\[CapitalSigma],a,\[Epsilon],\[Kappa],\[Tau]]/\[Delta]h2Ten[\[CapitalSigma],a,0,\[Kappa],\[Tau]]]]   
Ebend[\[CapitalSigma]_,\[Kappa]_,R_]:=(2\[Kappa])/R^2 \[CapitalSigma] ;
Eadh[\[CapitalSigma]_,\[Epsilon]_,a_]:=- \[CapitalSigma] \[Epsilon]/a^2;
Esurf[\[CapitalSigma]_,\[Tau]_,R_]:=+\[Tau] \[CapitalSigma]^2/(4\[Pi] R^2);
Sm[\[CapitalSigma]_,a_,\[Epsilon]_,\[Kappa]_,\[Tau]_]:=-(\[CapitalSigma]/a^2)  Log[Sqrt[\[Delta]h2Ten[\[CapitalSigma],a,\[Epsilon],\[Kappa],\[Tau]]/\[Delta]h2Ten[\[CapitalSigma],a,0,\[Kappa],\[Tau]]]]  ;


(* ::Input:: *)
(*a=1.145;*)
(*\[Epsilon]=0.7;*)
(*\[Tau]=0.001;*)
(*R=6.17; (* = (10+1)/2*2^(1/6), position of the adhesion energy minimum between nanoparticle and membrane bead *)*)


(* ::Input:: *)
(*Manipulate[Plot[{F[\[CapitalSigma],a,\[Epsilon],\[Kappa],\[Tau],R,\[Epsilon]]/.{\[Kappa]->10^\[Kappa]exp},(2\[Kappa])/R^2 \[CapitalSigma]/.{\[Kappa]->10^\[Kappa]exp} ,- \[CapitalSigma] \[Epsilon]/a^2,+\[Tau] \[CapitalSigma]^2/(4\[Pi] R^2),-(\[CapitalSigma]/a^2) kT Log[Sqrt[\[Delta]h2Ten[\[CapitalSigma],a,\[Epsilon],\[Kappa],\[Tau]]/\[Delta]h2Ten[\[CapitalSigma],a,0,\[Kappa],\[Tau]]]]/.{\[Kappa]->10^\[Kappa]exp}},{\[CapitalSigma],0,30},PlotLegends->{"Tot","\!\(\*SubscriptBox[\(E\), \(bend\)]\)","\!\(\*SubscriptBox[\(E\), \(adh\)]\)","\!\(\*SubscriptBox[\(E\), \(surf\)]\)","-TS"}, AxesLabel->{"\[CapitalSigma]"}],{\[Kappa]exp,-6,3}] *)
(* (* For a given \[Kappa] (adjustable through the sliding switch), plot energy contributions as a function of wrapping area \[CapitalSigma]*)*)


(* ::Input:: *)
(*\[Kappa]List=Catenate[{Table[10^\[Kappa]exp,{\[Kappa]exp,-3,0.5,0.1}],Table[10^\[Kappa]exp,{\[Kappa]exp,0.5,Log10[30],0.04}],Table[10^\[Kappa]exp,{\[Kappa]exp,Log10[30],Log10[40],0.01}],Table[10^\[Kappa]exp,{\[Kappa]exp,Log10[40],2,0.1}]}];   (* List of \[Kappa] points *)*)
(*Fminim=NMinimize[{F[\[CapitalSigma],a,\[Epsilon],#,\[Tau],R],\[CapitalSigma]>=0},\[CapitalSigma]]  & /@ \[Kappa]List//Chop; (* Minimize constrained free energy F(\[CapitalSigma]) *)*)
(*\[CapitalSigma]b\[Kappa]List=Table[{\[Kappa]List[[i]], \[CapitalSigma]/.Fminim[[i,2]]},{i,1,Length[\[Kappa]List]}];   (* Store list of Subscript[\[CapitalSigma], bound] that minimises free energy, for each \[Kappa] in \[Kappa]List *)*)
(*F\[Kappa]List=Table[{\[Kappa]List[[i]], Fminim[[i,1]]},{i,1,Length[\[Kappa]List]}];   (* Store list of equilibrium free energy F (resulting from minimisation with respect to \[CapitalSigma], for each \[Kappa] in \[Kappa]List *)*)


(* ::Input:: *)
(*ListLogLinearPlot[\[CapitalSigma]b\[Kappa]List,AxesLabel->{"\[Kappa]","\!\(\*SubscriptBox[\(\[CapitalSigma]\), \(bound\)]\)"},Joined->True]   (* Plot of wrapped area Subscript[\[CapitalSigma], bound] vs \[Kappa] *)*)
(*plotFLogLinear=ListLogLinearPlot[{F\[Kappa]List,{#[[1]],Ebend[#[[2]],#[[1]],R]}&/@\[CapitalSigma]b\[Kappa]List,{#[[1]],Eadh[#[[2]],\[Epsilon],a]}&/@\[CapitalSigma]b\[Kappa]List,{#[[1]],Esurf[#[[2]],\[Tau],R]}&/@\[CapitalSigma]b\[Kappa]List,{#[[1]],Sm[#[[2]],a,\[Epsilon],#[[1]],\[Tau]]}&/@\[CapitalSigma]b\[Kappa]List},AxesLabel->{"\[Kappa]",""},Joined->True,PlotLegends->{"F","\!\(\*SubscriptBox[\(E\), \(bend\)]\)","\!\(\*SubscriptBox[\(E\), \(adh\)]\)","\!\(\*SubscriptBox[\(E\), \(surf\)]\)","-TS"},PlotLabel->"\[Epsilon]="<>ToString[\[Epsilon]]<>", \[Tau]="<>ToString[\[Tau]]<>", R="<>ToString[R], PlotRange->All]    (* Plot of free energy (and its components) vs \[Kappa] *)*)
(*(* *)
(*Same plots in linear scale: *)
(**)
(*ListPlot[\[CapitalSigma]b\[Kappa]List,AxesLabel\[Rule]{"\[Kappa]","Subscript[\[CapitalSigma], bound]"},Joined\[Rule]True]*)
(*plotFLinear=ListPlot[{F\[Kappa]List,{#[[1]],Ebend[#[[2]],#[[1]],R]}&/@\[CapitalSigma]b\[Kappa]List,{#[[1]],Eadh[#[[2]],\[Epsilon],a]}&/@\[CapitalSigma]b\[Kappa]List,{#[[1]],Esurf[#[[2]],\[Tau],R]}&/@\[CapitalSigma]b\[Kappa]List,{#[[1]],Sm[#[[2]],a,\[Epsilon],#[[1]],\[Tau]]}&/@\[CapitalSigma]b\[Kappa]List},AxesLabel\[Rule]{"\[Kappa]",""},Joined\[Rule]True,PlotLegends\[Rule]{"F","Subscript[E, bend]","Subscript[E, adh]","Subscript[E, surf]","-TS"},PlotLabel\[Rule]"\[Epsilon]="<>ToString[\[Epsilon]]<>", \[Tau]="<>ToString[\[Tau]]<>", R="<>ToString[R], PlotRange\[Rule]Full]*)
(**)*)


(* ::Input:: *)
(*PrintList=Prepend[*)
(* Table[{*)
(*\[Kappa]List[[i]], *)
(*Fminim[[i,1]],*)
(*Ebend[\[CapitalSigma]b\[Kappa]List[[i,2]],\[Kappa]List[[i]],R],*)
(*Eadh[\[CapitalSigma]b\[Kappa]List[[i,2]],\[Epsilon],a], *)
(*Esurf[\[CapitalSigma]b\[Kappa]List[[i,2]],\[Tau],R], *)
(*Sm[\[CapitalSigma]b\[Kappa]List[[i,2]],a,\[Epsilon],\[Kappa]List[[i]],\[Tau]]*)
(*},{i,1,Length[\[Kappa]List]}]*)
(*,*)
(*{"#kappa", "F", "E_bending", "E_adhesion", "E_surface", "-TS"}*)
(*] // Chop*)
(**)
(*Export["DurotaxisGit_FreeEnergyContributions_eps"<>ToString[\[Epsilon]]<>"_tau"<>ToString[\[Tau]]<>"_R"<>ToString[R]<>"_a"<>ToString[a]<>".dat", PrintList]*)



