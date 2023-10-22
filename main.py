from loguru import logger
import pendulum
import os


def write_log( message:str, level:str = "INFO", save_file:bool = True, dir_path:str = "logs"):
    '''The `write_log` function logs a message with a specified level and saves it to a file if `save_file`
    is True.
    
    Parameters
    ----------
    message : str
        The `message` parameter is a string that represents the log message that you want to record. It can
        contain any information that you want to log, such as error messages, status updates, or debugging
        information.
        
    level : str
        The "level" parameter is used to specify the severity level of the log message. It can be one of
        the following values:
        - TRACE,
        - DEBUG,
        - INFO,
        - SUCCESS,
        - WARNING,
        - ERROR,
        - CRITICAL.\n
        
    save_file : bool, optional
        The `save_file` parameter is a boolean flag that determines whether the log message should be saved
        to a file or not. By default, it is set to `True`, meaning the log message will be saved to a file.
        
    dir_path : str, optional
        The `dir_path` parameter is a string that specifies the directory path where the log files will be
        saved. By default, it is set to "logs", which means that the log files will be saved in a directory
        named "logs" in the current working directory. However, you can provide a
    
    '''
    date = pendulum.now()
    file_name = f"{date.year}-{date.month}-{date.day}-{date.hour}-{date.minute}-{date.second}.log"
    
    if save_file:
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        logger.add(f"{dir_path}/{file_name}", format="{time} {level} {message}", level=level)
    
    logger.log(level, message)
    
