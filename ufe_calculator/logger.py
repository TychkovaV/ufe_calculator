from enum import IntEnum


class LogLevel(IntEnum):
    INFO = 1
    WARNING = 2
    ERROR = 3
    def __str__(self) -> str:
        return {LogLevel.INFO : "INFO",
                LogLevel.WARNING : "WARNING",
                LogLevel.ERROR : "ERROR"}[self]


class Logger:
    def __init__(self, logfile="", clear=True, filtering_level=LogLevel.INFO):
        self.filter:LogLevel = filtering_level
        self.output:str = logfile
        if clear:
            self.clear()

    def clear(self) -> None:
        if self.output!="":  
            tmp = open(self.output, "rt")
            tmp.close() 

    def log(self, text:str, level:LogLevel) -> None:
        if level<self.filter:
            return
        if self.output=="":
            print(str(level)+": "+text)
        else:
            with open(self.output, "at") as f:
                f.write(str(level)+": "+text)
    def info(self, text:str) -> None:
        self.log(text=text, level=LogLevel.INFO)
    def warning(self, text:str) -> None:
        self.log(text=text, level=LogLevel.WARNING)
    def error(self, text:str) -> None:
        self.log(text=text, level=LogLevel.ERROR)    