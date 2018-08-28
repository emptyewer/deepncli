

def elapsed_time(start, stop):
    temp = stop - start
    hours = temp // 3600
    temp = temp - 3600 * hours
    minutes = temp // 60
    seconds = temp - 60 * minutes
    return hours, minutes, seconds
