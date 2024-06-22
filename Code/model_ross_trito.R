library(deSolve)
library(reshape2)

###One path, no dispersal- simplest form.


model_ross_trito <-  function(t, state, param) {
        with(as.list(c(state, param)),{
                
                ###Demographic parameters
                b_H <- param["b_H"] #Human birth rate
                b_P <- param["b_P"] #P.vector birth rate
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
                dH_S <- b_H * (N_H) - (FOI_P * H_S * (P_I/N_P)) - (FOI_S *H_S*(S_I/N_S))- (mu_H*H_S)
                
                ###Infected
                dH_I <- (FOI_P * H_S * (P_I/N_P)) + (FOI_S *H_S*(S_I/N_S))- (gamma*H_I) -(mu_H*H_I)
                
                ###Recovered
                dH_R <- (gamma*H_I)-(mu_H*H_R)
                
                
                ###P. vector
                ###Susceptible
                dP_S <-  (b_P* N_P) - (FOI_H_P *P_S* (H_S/N_H)) - mu_P *P_S - c_SP*(P_S)*(N_S)
                
                ###Infected
                dP_I <- (FOI_H_P * P_S * (H_S/N_H)) - mu_P*P_I- c_SP*(P_I)*(N_S)
                 
                
                ###S. vector 
                ###Susceptible
                dS_S <-  (b_S * N_S) - (FOI_H_S *S_S* (H_S/N_H)) - mu_S *S_S -  c_PS*(S_S)*(N_P)
                
                ###Infected
                dS_I <- (FOI_H_S * S_S * (H_S/N_H)) - mu_S*S_I -  c_PS*(S_I)*(N_P)
                
                return(list(c(dH_S,dH_I,dH_R,dP_S, dP_I, dS_S, dS_I)))
        }
        )
}

###BAD PARAM VALUES 
parameters_n <- c(
        b_H = 0.005, #Human birth rate
        b_P = 1, #P.vector birth rate
        b_S = 1, #S. vector birth rate
        mu_H = 0.005, #Human death rate
        mu_P = 1, #P. vector death rate
        mu_S = 1, #S. vector death rate
        
        a_P =3, #biting rate of the p. vector
        a_S = 3, #biting rate of the s.vector
        
        phi_P = 0.90, #transmission probability of p. vector
        phi_S = 0.10, #transmission probability of s. vector
        phi_H  = 0.90, #transmission probability of human
        
        # Recovery rate
        gamma = 1/7,  #recovery rate of infected human
        
        #competition coefficient
        c_PS = 0.001, #competitition effect of p.vector on s.vector
        c_SP = 0  #competitition effect of s.vector on p.vector
        )

inits_n <- c(H_S = 900,H_I =100 ,
             H_R = 0 ,P_S = 900, P_I= 100, S_S = 900, S_I= 100)

times <- seq(0,25, by = 0.5)


out_DDE <- data.frame(ode(y = inits_n, times = times, func = model_ross_trito,
                          parms = parameters_n, method = "lsoda"))


melted_out_DDE <- melt(out_DDE, id.vars = 'time')

melted_out_DDE$grouping <- NA
melted_out_DDE $grouping<- ifelse(melted_out_DDE$variable %in% c("H_S", "H_I", "H_R"),
                                  "H",  melted_out_DDE $grouping)
melted_out_DDE $grouping<- ifelse(melted_out_DDE$variable %in% c("P_S", "P_I"),
                                  "P",  melted_out_DDE $grouping)
melted_out_DDE $grouping<- ifelse(melted_out_DDE$variable %in% c("S_S", "S_I"),
                                  "S",  melted_out_DDE $grouping)


ggplot(melted_out_DDE , aes(x= time, y= log10(value+1), color = variable))+
        geom_line(size = 1)+facet_wrap(~grouping)+
        scale_color_viridis(discrete = TRUE, option = 'turbo')+
        theme_classic()+
        xlab("Time")+
        ylab("Log10 abundance")+
        theme(strip.background = element_blank())
        

