def hms_string(sec_elapsed) -> str:
    """Nicely formatted time string.
    :param sec_elapsed time in ms
    :return
        hh:mm:ss.msms
    """
    h = int(sec_elapsed / (60 * 60))
    m = int((sec_elapsed % (60 * 60)) / 60)
    s = sec_elapsed % 60
    return "{}:{:>02}:{:>05.2f}".format(h, m, s)
