#!/bin/bash

# --------------------------------------------------------------------------------
# A set of functions to post to gomez-chiarri-lab slack

# Refer to - 
# - https://gist.github.com/dopiaza/6449505
# - https://api.slack.com/incoming-webhooks

# this URL is only valid to post in #cluster-jobs channel
CLUSTER_JOBS_SLACK_URL="https://hooks.slack.com/services/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"

FORMATTER_SCRIPT="$(dirname ${BASH_SOURCE[0]})/slack_message_formatter.py"

function post_slack_message()  {
    # two inputs - <channel name> and <message to be posted>
    # e.g.
    #   post_slack_message cluster-jobs "hi there"
    
    case $1 in
    "cluster-jobs")
        url="${CLUSTER_JOBS_SLACK_URL}"
        ;;
    *)
        echo "ERROR: Unsupported slack channel '$1'"
        return 1
        ;;
    esac

    if [ -z "$2" ]; then
        echo "ERROR: No message given"
        return 1
    fi

    local from="cluster-job-bot"
    if [ "$3" ]; then
        from=$3
    fi

    local message_file_name=/tmp/slack.msg.$$
    if [ "$2" == "-" ] || [ "${2:0:1}" == "@" ]; then
        # message is a file
        cat "$2" | python $FORMATTER_SCRIPT $from > ${message_file_name}
    else
        echo "$2" | python $FORMATTER_SCRIPT $from > ${message_file_name}
    fi

    # variable $? indicates whether the last command succeed (i.e. equals 0) or not
    # only on successfully writing the message to $message_file_name in json format the
    # way slack wants, should we use curl to send it to slack
    if [ "$?" == 0 ]; then
        echo curl -X POST -H 'Content-type: application/json' --data "@${message_file_name}" $url
        curl -X POST -H 'Content-type: application/json' --data "@${message_file_name}" $url
    fi

    # if succeeded this function also deemed success - returning 0
    return $?
}


