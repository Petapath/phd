def voltage_lin_sweep(Vi,s,t):
    return (Vi + s*t)

def voltage_bidirectional_sweep(Vi,s,switch_time,t):
    if(t>=switch_time):
        ret = (Vi + s*switch_time) - s*(t-switch_time)
    else:
        ret = (Vi + s*t)
    return ret

# this function, when used in the context of odeint 
# is not a vector function, t is simply a scalar
#
def voltage_pulse(Vmin,Vmax,fr,to,t):
    if( t >= fr and t<to ):
        return Vmax
    else:
        return Vmin
