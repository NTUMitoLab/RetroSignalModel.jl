var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = RetroSignalModel","category":"page"},{"location":"#RetroSignalModel","page":"Home","title":"RetroSignalModel","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RetroSignalModel]","category":"page"},{"location":"#RetroSignalModel.actM1","page":"Home","title":"RetroSignalModel.actM1","text":"Activation Model:\n\n\n\n\n\n","category":"constant"},{"location":"#RetroSignalModel.rtgM1","page":"Home","title":"RetroSignalModel.rtgM1","text":"Rtg Model: version 1\n\n\n\n\n\n","category":"constant"},{"location":"#RetroSignalModel.rtgM2","page":"Home","title":"RetroSignalModel.rtgM2","text":"This model is modified version of rtgM1. Major changes:\n\nChange input to be michaelis menten\nChange the translocation to be Michaelis menten\nRtg1p and Rtg3p enter nucleus separately\n\n\n\n\n\n","category":"constant"},{"location":"#RetroSignalModel.rtgM3","page":"Home","title":"RetroSignalModel.rtgM3","text":"This model is modified version of rtgM2. Major changes:\n\nChange input to be michaelis menten\nChange the translocation to be Michaelis menten: (which is not ture). See the supplemental material of https://www.pnas.org/content/pnas/suppl/2006/06/30/0604085103.DC1/04085SuppText.pdf\nRtg1p and Rtg3p enter nucleus separately\n\nRef:\n\nhttps://science.sciencemag.org/content/339/6118/460\nhttps://www.pnas.org/content/pnas/suppl/2006/06/30/0604085103.DC1/04085SuppText.pdf\nhttps://www.cell.com/molecular-cell/pdf/S1097-2765(04)00179-0.pdf (Reaction of Bmh-Mks on Rtg13 dimer)\n\n\n\n\n\n","category":"constant"},{"location":"#RetroSignalModel.rtgM4","page":"Home","title":"RetroSignalModel.rtgM4","text":"Modified version from rtgM3\n\nChange input to Hill kinetics\n\n\n\n\n\n","category":"constant"},{"location":"#RetroSignalModel.solve_dym_rtgM4-NTuple{4,Any}","page":"Home","title":"RetroSignalModel.solve_dym_rtgM4","text":"ODE function of rtgM4  Created and modified via ModelConv.createfuncexp(model)\n\nArgumemt\n\n'u::Vector': systemic component of model\n'p::Vector': parameters\n'tspan::Tuple(1,2)': Start time , and ending\n'sig::Function': sig(t) is the signal generator. The signal should be differentiable.\n\nExamples\n\n> model = Model.rtgM4\n> solution = get_solution(;i=4)\n> u,p = rtgPar()\n> sig = Waves.SIN(amp=3, freq=0.5/(2*pi), phi= 2*pi * 3/4)\n> tspan = (0.0,200.0)\n> sol = solve_ode(param.u, param.p, tspan, sig)\n\nNote\n\nIf the parameter set is not in steady-state when t=0 with sig(t). This function automatically calculate the steady-state and start the simulation. However, it takes seconds to\n\n\n\n\n\n","category":"method"}]
}
