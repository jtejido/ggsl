package err

import (
	"log"
)

type (
	ErrorHandlerType  = func(reason, file string, line, gsl_errno int)
	StreamHandlerType = func(label, file string, line int, reason string)
)

var (
	errorHandler ErrorHandlerType = nil
)

type GSLError interface {
	Status() int
	Error() string // should embed 'error' interface, but running into errors compiling on Windows machines.
}

type gslerror struct {
	status  int
	message string
}

func New(status int, text string) GSLError {
	return &gslerror{status, text}
}

func (err *gslerror) Status() int {
	return err.status
}

func (err *gslerror) Error() string {
	return err.message
}

func Error(reason, file string, line, gsl_errno int) {
	if errorHandler != nil {
		errorHandler(reason, file, line, gsl_errno)
		return
	}

	StreamPrintf("ERROR", file, line, reason)
	log.Printf("Default GGSL error handler invoked.\n")
	panic(reason)
}

func SetErrorHandler(new_handler ErrorHandlerType) ErrorHandlerType {
	previous_handler := errorHandler
	errorHandler = new_handler
	return previous_handler
}

func SetErrorHandlerOff() ErrorHandlerType {
	previous_handler := errorHandler
	errorHandler = NoErrorHandler
	return previous_handler
}

func NoErrorHandler(reason, file string, line int, gsl_errno int) {
	/* do nothing */
	return
}
