#!/usr/bin/env python
#-*- coding: utf-8 -*-

import itchat
@itchat.msg_register(itchat.content.TEXT)
def print_content(msg):
    return msg.user.send('【自动回复】对不起，您联系的用户暂时不在微信旁。如有急事，请联系：15600927092；jerryshen68@hotmail.com')
itchat.auto_login(hotReload=True)#,enableCmdQR = True)
itchat.run()
























