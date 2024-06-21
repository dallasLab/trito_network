model_ross_trito <-  function(t, state, param) {
        with(as.list(c(state, param)),{
                
                ###Demographic parameters
                b_H <- param["b_H"] #Human birth rate
                b_P <- param["b_V"] #P.vector birth rate
                b_S <- param["b_S"] #S. vector birth rate
                mu_H <- param["mu_H"] #Human death rate
                mu_P <- param["mu_P"] #P. vector death rate
                mu_S <- param["mu_S"] #S. vector death rate
                
                #Force of infection parameters
                a_P <- param["a_P"] #biting rate of the p. vector
                a_S <- param["a_S"] #biting rate of the s.vector
                
                phi_P <- param["phi_P"] #transmission probability of p. vector
                phi_S <- param["phi_S"] #transmission probability of s. vector
                phi_H  <- param["phi_H"] #transmission probability of human
                
                # Recovery rate
                gamma <- param["gamma"] #recovery rate of infected human
                
                #competition coefficient
                c_PS <- param["c_PS"] #competitition effect of p.vector on s.vector
                c_SP <- param["c_SP"] #competition effect of s.vector on p.vector
                
                ### FOI
                FOI_P <- a_P * phi_P #FOI for a primary vector
                FOI_S <- a_S * phi_S #FOI for a secondary vector
                FOI_H_P <- a_P * phi_H #FOI for a human to a primary vector
                FOI_H_S <- a_S * phi_S #FOI for a human to a secondary vector
                
                ### Population size
                N_H <- H_S + H_I + H_R #human host population 
                N_P <- P_S + P_I #p. vector population
                N_S <- S_S + S_I #s. vector population
                
                ###Human host 
                
                ###Susceptible 
                dH_S <- b_H - (FOI_P * H_S * (P_I/N_P)) - (FOI *H_S*(S_I/N_S))- (mu_H*H_S)
                
                ###Infected
                dH_I <- (FOI_P * H_S * (P_I/N_P)) + (FOI *H_S*(S_I/N_S))- (gamma*H_I) -(mu_H*H_I)
                
                ###Recovered
                dH_R <- (gamma*H_I)-(mu_H*H_R)
                
                
                ###P. vector
                ###Susceptible
                dP_S <-  b_P - (FOI_H_P *P_S* (H_S/N_H)) - mu_P *P_S
                
                ###Infected
                dP_I <- (FOI_H_P * P_S * (H_S/N_H)) - mu_P*P_I
                 
                
                ###S. vector 
                ###Susceptible
                dS_S <-  b_S - (FOI_H_S *S_S* (H_S/N_H)) - mu_S *S_S
                
                ###Infected
                dS_I <- (FOI_H_S * S_S * (H_S/N_H)) - mu_S*S_I
                
                return(list(c(dH_S,dH_I,dH_R,dP_S, dP_I, dS_S, dS_I)))
        }
        )
}


###BAD PARAM VALUES 
parameters_n <- c(
        b = 100,
        alphaH =1/10,
        beta = 0.005,
        delta = 0.1,
        alphaP = 1/10,
        gamma = 0,
        nH = 3,
        nP = 3)

inits_n <- c(HJ = rep(0,parameters_n['nH']), 
             HA = 0,
             PJ =  rep(0,parameters_n['nP']),
             PA  = 0)

times <- seq(0, 1000, by = 1)


out_DDE <- data.frame(ode(y = inits_n, times = times, func = PP_variance ,
                          parms = parameters_n, method = "lsodar"))


