#!/usr/bin/env python
#-*- coding: utf-8 -*-

import itchat
to_place=input('Where to go: ')
if to_place in ['健身房','0'] :
    send_message='【自动回复】对不起，我去健身了，暂时无法回复。'
else:
    send_message='【自动回复】对不起，您联系的用户暂时不在微信旁。如有急事，请联系：15600927092；jerryshen68@hotmail.com'
@itchat.msg_register(itchat.content.TEXT)
def print_content(msg):
    return msg.user.send(send_message)
itchat.auto_login(hotReload=True)#,enableCmdQR = True)
itchat.run()
























