//bayes.hh
/// General structures for Bayesian data analysis
/// At initial writing the plan for this is in flux. Some elements we
/// probably want here are:
///  -Sampler:     A generalized class for a Bayesian sampler, which could be track to various MCMC codes, or a 
///                Multinest code, for instance.
///  -DataModel:   A generalized class for modeling related to an instrument, this would generally include everything
///                needed for computing a likelihood that aligns with the data.  This would include the instrumental
///                model, and the stochasitic model for the likelihood. Some of this, such as the basics of chisq
///                computation, are generic and should belong in the base class, or close derivatives.  Other methods
///                such as how to read data files or instrumental parameters go in concrete derived classes. 
///  -SignalModel: Here we include specifics of the signal model including incident data simulation
///                
