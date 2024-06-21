###Modeling metapopulation 

model_ross_trito_metapopulation <- function(t, state, param, patch_no, disp.contact) {
        
        with(as.list(c(state, param)), {
                
        num_patch = patch_no
        
        ###The metapopulation
        H_S = matrix(state[1:num_patch], ncol = 1)
        H_I = matrix(state[(num_patch+1):(2*num_patch)], ncol = 1)
        H_R = matrix(state[((2*num_patch)+1):(3*num_patch)], ncol = 1)
        P_S = matrix(state[((3*num_patch)+1):(4*num_patch)], ncol = 1)
        P_I = matrix(state[((4*num_patch)+1):(5*num_patch)], ncol = 1)
        S_S = matrix(state[((5*num_patch)+1):(6*num_patch)], ncol = 1)
        S_I = matrix(state[((6*num_patch)+1):(7*num_patch)], ncol = 1)
        
        ###Demographic parameters
        b_H <- matrix(rep(param["b_H"],num_patch) , ncol = 1) #Human birth rate
        b_P <- matrix(rep(param["b_P"],num_patch) , ncol = 1)  #P.vector birth rate
        b_S <-  matrix(rep(param["b_S"],num_patch) , ncol = 1)  #S. vector birth rate
        mu_H <- matrix(rep(param["mu_H"],num_patch) , ncol = 1)  #Human death rate
        mu_P <-  matrix(rep(param["mu_P"],num_patch) , ncol = 1)  #P. vector death rate
        mu_S <- matrix(rep(param["mu_S"],num_patch) , ncol = 1)  #S. vector death rate
        
        #Force of infection parameters
        a_P <- matrix(rep(param["a_P"],num_patch) , ncol = 1)  #biting rate of the p. vector
        a_S <-  matrix(rep(param["a_S"],num_patch) , ncol = 1)  #biting rate of the s.vector
        
        phi_P <-  matrix(rep(param["phi_P"],num_patch) , ncol = 1)  #transmission probability of p. vector
        phi_S <-  matrix(rep(param["phi_S"],num_patch) , ncol = 1)   #transmission probability of s. vector
        phi_H  <-  matrix(rep(param["phi_H"],num_patch) , ncol = 1)   #transmission probability of human
        
        # Recovery rate
        gamma <-  matrix(rep(param["gamma"],num_patch) , ncol = 1)   #recovery rate of infected human
        
        #competition coefficient
        c_PS <-  matrix(rep(param["c_PS"],num_patch) , ncol = 1)   #competitition effect of p.vector on s.vector
        c_SP <-  matrix(rep(param["c_SP"],num_patch) , ncol = 1)  #competition effect of s.vector on p.vector
        
        ### FOI
        FOI_P <- matrix(a_P * phi_P, ncol = 1) #FOI for a primary vector
        FOI_S <- matrix(a_S * phi_S , ncol = 1) #FOI for a secondary vector
        FOI_H_P <-  matrix(a_P * phi_H,ncol =1) #FOI for a human to a primary vector
        FOI_H_S <- matrix(a_S * phi_S,ncol = 1) #FOI for a human to a secondary vector
        
        
    
        
        ### Population size
       N_H <- matrix(H_S + H_I + H_R)#human host population 
       N_P <- matrix(P_S + P_I) #p. vector population
       N_S <- matrix(S_S + S_I) #s. vector population
        
        ###Human host 
        ###Susceptible 
        dH_S <- b_H * (N_H) - (FOI_P * H_S * (P_I/N_P)) - (FOI_S *H_S*(S_I/N_S))- (mu_H*H_S) 
        

        ###Infected
        dH_I <- (FOI_P * H_S * (P_I/N_P)) + (FOI_S *H_S*(S_I/N_S))- (gamma*H_I) -(mu_H*H_I)
        
        ###Recovered
        dH_R <- (gamma*H_I)-(mu_H*H_R)
        
        ###There are no dispersal here
        
        
        ###P. vector
        ###Susceptible
        dP_S <-  (b_P* N_P) - (FOI_H_P *P_S* (H_S/N_H)) - mu_P *P_S - c_SP*(P_S)*(N_S) +
                (t(disp.contact) %*% P_S) -  (disp.contact %*% P_S)
        
        ###Infected
        dP_I <- (FOI_H_P * P_S * (H_S/N_H)) - mu_P*P_I- c_SP*(P_I)*(N_S) + 
                (t(disp.contact) %*% P_I) -  (disp.contact %*% P_I)
        
        
        ###S. vector 
        ###Susceptible
        dS_S <-  (b_S * N_S) - (FOI_H_S *S_S* (H_S/N_H)) - mu_S *S_S -  c_PS*(S_S)*(N_P)+
                (t(disp.contact) %*% S_S) -  (disp.contact %*% S_S)
        
        
        ###Infected
        dS_I <- (FOI_H_S * S_S * (H_S/N_H)) - mu_S*S_I -  c_PS*(S_I)*(N_P) + 
                (t(disp.contact) %*% S_I) -  (disp.contact %*% S_I)
        
        return(list(c(dH_S,dH_I,dH_R,dP_S, dP_I, dS_S, dS_I)))
        }
        )
}

###BAD PARAM VALUES 
parameters_n <- c(
        b_H = 0.005, #Human birth rate
        b_P = 100, #P.vector birth rate
        b_S = 100, #S. vector birth rate
        mu_H = 0.005, #Human death rate
        mu_P = 100, #P. vector death rate
        mu_S = 100, #S. vector death rate
        
        a_P =3, #biting rate of the p. vector
        a_S = 3, #biting rate of the s.vector
        
        phi_P = 0.90, #transmission probability of p. vector
        phi_S = 0.10, #transmission probability of s. vector
        phi_H  = 0.40, #transmission probability of human
        
        # Recovery rate
        gamma = 1/7,  #recovery rate of infected human
        
        #competition coefficient
        c_PS = 0, #competitition effect of p.vector on s.vector
        c_SP = 0  #competitition effect of s.vector on p.vector
)

inits_n <- c(H_S = rep(900,3), H_I =rep(100 ,3),
             H_R = rep(0,3) ,
             P_S = rep(15000,3), 
             P_I= rep(0,3),
             S_S = rep(15000,3),
             S_I= rep(0,3))

df.contact <- matrix(c(0,1,1,
                     1,0,0,
                     0,1,0), ncol =3)
        

out <- ode(y=inits_n , times=times, patch_no=3, func=model_ross_trito_metapopulation,
           parms=parameters_n,disp.contact=df.contact)

colnames(out)
plot(out[,'time'],out[,14],type= 'l')
lines(out[,'time'],out[,15],col = 'red')
lines(out[,'time'],out[,16],col = 'green')
lines(out[,'time'],out[,20],col = 'blue',type='l')
lines(out[,'time'],out[,21],col = 'orange', size = 2)
