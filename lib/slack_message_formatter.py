import sys
import json

text = sys.stdin.read()
if len(sys.argv) > 1:
    user = sys.argv[1]
else:
    user = "bot"

msg = { 'mrkdwn' : True, 'username' : user,
        'text' : "```" + text + "```" }
json.dump(msg, sys.stdout)
