package err

import (
	"log"
	"os"
)

var (
	stream        *os.File          = nil
	streamHandler StreamHandlerType = nil
)

func StreamPrintf(label, file string, line int, reason string) {
	if stream == nil {
		stream = os.Stderr
		log.SetOutput(stream)
	}

	if streamHandler != nil {
		streamHandler(label, file, line, reason)
		return
	}

	log.Printf("ggsl: %s:%d: %s: %s\n", file, line, label, reason)

}

func SetStreamHandler(new_handler StreamHandlerType) StreamHandlerType {
	previous_handler := streamHandler
	streamHandler = new_handler
	return previous_handler
}

func SetStream(new_stream *os.File) *os.File {
	if stream == nil {
		stream = os.Stderr
	}

	previous_stream := stream
	stream = new_stream
	return previous_stream
}
