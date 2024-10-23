import numpy as np

if __name__ == '__main__':

    weather_file = '/home/johnson.7080/gulls_sj/weather/WFIRSTLSST1-72.weather'
    init_zero = 3.
    final_zero = 3.
    length = 72.

    total_len = init_zero+length+final_zero
    
    days = np.linspace(0,total_len,total_len*4+1,endpoint=True)


    with open(weather_file,'w') as w_file:
        for day in days:
            if day>=init_zero and day<(total_len-final_zero):
                w_file.write("%0.2f %i\n"%(day,1))
            else:
                w_file.write("%0.2f %i\n"%(day,0))

    exit()
        
